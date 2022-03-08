rm(list=ls())
library(tidyverse)
library(data.table)
library(readtext)
library(tidytext)
library(ngram)

dir.input <- "path/to/notes"
options(datatable.fread.datatable=F)
true.annot <-
  fread(file.path(dir.input,"annotations.csv")) %>%
  select(starts_with(c("file","caregiver presence","home","inst","formal","informal"))) %>%
  mutate(across(starts_with(c("caregiver presence","home","inst","formal","informal")), as.numeric))
names(true.annot) <- c("doc_id","presence","home","institution","formal","informal")
true.annot <- true.annot %>% arrange(doc_id)
files.train <-
  list.files(dir.input, full.names=T, pattern=".txt")
data.train <- readtext(files.train) %>%
  distinct(doc_id, .keep_all=T) %>%
  arrange(doc_id) %>%
  mutate(doc_id=str_remove(doc_id, ".txt"))
dict <- fread("path/to/dictionary.csv") %>%
  distinct(ngram, .keep_all=T) %>% arrange(ngram)
names(dict) <- c("ngram","presence","home","institution","formal","informal")
names(true.annot) <- c("doc_id",
                       dict %>% select(-contains("gram")) %>% names)
true.annot <- true.annot %>% distinct(doc_id, .keep_all=T) %>%
  arrange(doc_id) %>% mutate(doc_id=str_remove(doc_id, ".txt"))

# Set parameters
ngram_count <- 6 #Up to how many ngrams should be tokenized?
discharge_tol <- 3 #How many words following or preceding a dictionary term should be checked for "discharge"
patient_tol <- 4 #How many words following or preceding a dictionary term should be checked for 'patient' logic
patient_verbs <- unname(unlist(fread("path/to/patient_verbs.csv")))
# Tokenization
ngramize <- function(data, n) {
  lapply(1:n, function(i)
    data %>% unnest_ngrams(ngram, text, n=i)) %>% bind_rows
}

#Iterate through each unique medical note
preds <- lapply(data.train %>% pull(doc_id), function(id) {
  pred <- data.frame(ngram=character(),
                     presence=numeric(),
                     home=numeric(),
                     institution=numeric(),
                     formal=numeric(),
                     informal=numeric())
  data.tmp <- data.train %>% filter(doc_id==id) %>%
    mutate(text=str_replace_all(text, "\\bPT\\b", "physical therapy")) %>%
    mutate(text=str_to_lower(text))
  # remove/replace select patterns in notes
  patterns.rm <-
    c("nurse navigator","navigator nurse",
      "patient portal", "patient name",
      "relationship to patient:?\\s*\\n",
      "verbalize(s|d) understanding",
      "e\\.?g\\.?,?\\s?visiting nurses?",
      "patient(\\s?&\\s?| or | &/or )caregiver",
      "family (medicine|practice|doctor|physician)",
      "alone with family",
      "family history.*\\n") # "FAMILY HISTORY:mother"
  patterns.rp <- c("want(s?|ing|ed) [[:alpha:]]* (patient|pt)",
                   "want(s?|ing|ed) (patient|pt)")
  data.tmp <- data.tmp %>%
    mutate(text=str_remove_all(text,paste(patterns.rm,collapse="|")),
           text=str_replace_all(text,paste(patterns.rp,collapse="|"),"patient"),
           text=str_replace_all(text,"pt or ot|pt and ot|pt/ot",
                                "physical therapy and ot"),
           text=str_replace_all(text,"e-mail","email"),
           text=str_replace_all(text,"patient(s|'s|\\ns|\\ts|\\ss|)\\spartner",
                                "patient partner"))
  # look for "visit" in the neighborhood of "nurse" or "nurses"
  VN <- data.tmp %>% pull(text) %>%
    str_detect("visit(ing|)[^.:\n\t]*nurses?|nurses?[^.:\n\t]*\\.?[^.:\n\t]*visit")
  if (VN) pred <- pred %>% add_row(dict %>% filter(ngram=="visiting nurse"))
  # drop sleeping/hearing aid(s)
  if (grepl("(sleep|sleeping|hear|hearing) aids?", data.tmp %>% pull(text)))
    data.tmp <- data.tmp %>%
    mutate(text=str_remove_all(text, "\\baids?"))
  
  ngrams <- data.tmp %>% ngramize(ngram_count)
  ngrams.inst <- intersect(ngrams %>% pull(ngram),
                           dict %>% filter(institution==1) %>% pull(ngram))
  
  sents <- data.tmp %>% pull(text) %>% strsplit("!|\\.|\\?|\\n+") %>% unlist
  # drop terms in sentences with institutional ngrams
  if (length(ngrams.inst)>0) {
    terms.drop <-
      c("patient","pt","care giver","caregiver","guardian")
    include.inst.ngram <-
      sapply(sents, function(sent) grepl(paste(ngrams.inst, collapse="|"), sent))
    invisible(sapply(which(include.inst.ngram),
                     function(id) sents[id] <<-
                       str_remove_all(sents[id],paste(terms.drop, collapse="|"))))
    data.tmp <- data.tmp %>%
      mutate(text=paste(sents, collapse=""))
  }
  # drop institutional ngrams as evidence for potential/past/declined caregivers
  patterns.inst <- c("return(ed|ing|)[^.]*from.*",
                     "discharg(ed|ing|es?)[^.]*from.*",
                     "[^.\n\t]*wait(ing\\s|\\s|)list[^.\n\t]*",
                     "[^.\n\t]*\\bcancel[^.\n\t]*",
                     "[^.\n\t]*\\bdeclin[^.\n\t]*",
                     "[^.\n\t]*\\brehab[^.\n\t]*",
                     "[^.\n\t]*\\bapprov[^.\n\t]*",
                     "[^.\n\t]*\\brequir[^.\n\t]*",
                     "[^.\n\t]*\\bsuggest[^.\n\t]*",
                     "[^.\n\t]*\\brequest[^.\n\t]*")
  sents.inst.false <- sents %>% str_extract(paste(patterns.inst,collapse="|"))
  sents.inst.false <- sents.inst.false[!is.na(sents.inst.false)]
  if (length(sents.inst.false)>0 & length(ngrams.inst)>0) {
    ngrams.inst.kept <- ngrams.inst[!sapply(ngrams.inst, function(ngram)
      sents.inst.false %>% str_detect(ngram) %>% any)]
    if (length(ngrams.inst.kept)>0)
      pred <- pred %>% add_row(dict %>% filter(ngram %in% ngrams.inst.kept))
  }
  ngrams <- data.tmp %>% ngramize(ngram_count)
  # look for evidence of self-care
  onegrams <- data.tmp %>% unnest_ngrams(ngram, text, n=1)
  self.care <-
    sapply(onegrams %>% pull(ngram) %in% c("patient","pt","patient's") %>% which,
           function(idx) {
             neighborhood <- onegrams %>%
               slice(max(idx-patient_tol,1):min(idx+patient_tol,NROW(onegrams))) %>%
               pull(ngram)
             sapply(patient_verbs, grepl, x=neighborhood) %>% any
           }) %>% any
  if (data.tmp %>%
      str_detect("relationship to patient\\s?:?\\s*(self|patient|pt)"))
    self.care <- T
  if (self.care) {
    pred <- pred %>%
      add_row(setNames(data.frame("patient",0,1,0,0,0),
                       c("ngram","presence","home","institution","formal","informal")))
  }
  # non-institutional ngrams
  ngrams.noninst <- intersect(ngrams %>% pull(ngram),
                              dict %>% filter(institution!=1) %>% pull(ngram))
  if (length(ngrams.noninst)>0) {
    pred <- pred %>% add_row(dict %>% filter(ngram %in% ngrams.noninst))
    if (exists("ngrams.inst.kept")) {
      ngrams.inst.kept <- intersect(ngrams.inst.kept,
                                    dict %>% filter(institution!=1) %>% pull(ngram))
      if (length(ngrams.inst.kept)>0) { # drop the pick-ups of home=1
        pred <- pred %>%
          filter(!ngram%in%c("patient"))
        ngrams.inst.kept <- NULL
      }
    }
  }
  # summarize
  if (NROW(pred)>0) {
    res <- pred %>% select(ngram) %>%
      summarise(ngram=paste(ngram, collapse=", ")) %>%
      bind_cols(pred %>% select(-ngram) %>% summarise(across(.fns=max)))
    # if institution is 1, home should be 0
    if (res %>% pull(institution)==1)
      res <- res %>% mutate(home=0)
  } else {
    res <- setNames(data.frame(NA,0,0,0,0,0),
                    c("ngram","presence","home","institution","formal","informal"))
  }
  return(res %>% mutate(doc_id=id) %>% relocate(doc_id))
}) %>% bind_rows