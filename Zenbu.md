Upload to the Zenbu genome browser.
===================================

This page shows how we loaded the aligned data on the "[Zenbu][]" genome
browser.  It is provided as a simple Markdown file for reference because the
upload had to be done only once.  We are currently improving and generalising
the commands below ease and automate browser uploads with knitr.

[Zenbu]: http://fantom.gsc.riken.jp/zenbu

For these commands to work, the user needs to have configured the tool
`zenbu_upload` (available from Zenbu version 2.9.1).  Here we upload on RIKEN's
main public server, in a collaboration that we use to privately manage our open
data.  In order to reuse the code below, one would need to create one own's
collaboration, or to directly upload to the "public" collaboration.

```{r}
ZENBU_COLLAB <- "BLFNw_m6NRVgdC2XaT2NcB"
```

Accessory functions
-------------------

Ad-hoc wrapper to the shell command `zenbu_upload`

```{r}
zenbuUpload <- function ( ...
                        , URL="http://fantom.gsc.riken.jp/zenbu"
                        , verbose=FALSE
                        , echo=FALSE
                        , stdout=TRUE) {
  zenbu <- 'zenbu_upload'
  url <- c('-url', URL)
  args <- sapply(c(url, ...), shQuote)
  if (verbose == TRUE) print(paste(c(zenbu, args), collapse=' '))
  if (echo    == FALSE) {
      system2(zenbu, args, stdout=stdout)
  } else {
      system2('echo', c(zenbu, args), stdout=stdout)
  }
}
```

To convert from file name to sample name.

```{r}
bedToSample <- function(BED)
  BED %>%
    sub("RunA", RunA, .) %>%
    sub("RunB", RunB, .) %>%
    sub(".bed", "", .) %>%
    sub("_R1", "", .)
```

To produce a string of keywords that will uniquely identify a sample.

```{r}
standardDescription <- function(BED)
  paste( LIBRARY
       , BED %>% bedToSample
       , "C1-CAGE-preview"
       , "CAGEscan_fragments"
       , "KnitrUpload")
```

To upload CAGEscan fragments, only if they have not yet been uploaded.

```{r}
safeZupload <- function (BED) {
  # Warning: the "notFound" function ignores queued uploads...
  notFound <- function(BED)
    zenbuUpload( "-list", "-filter", standardDescription(BED)) %>%
      tail(1) %>%
      grepl ("0 uploads --- 0 featuresources --- 0 experiments --- \\[0 total sources\\]", .)
  if(notFound(BED))
    zenbuUpload( "-file",        paste0("output/cagescan_fragments/", BED)
               , "-name",        BED  %>% sub(".bed", "", .)
               , "-assembly",    "hg38"
               , "-desc",        standardDescription(BED)
               , "-collab_uuid", ZENBU_COLLAB
               , "-singletag_exp"
               , stdout="")
}
```

Functions to add metadata tags.

```{r}
zenbuTag <- function (filter, key, value, mode='add', ...) {
  args <- c('-mdedit', filter, mode, key, value)
  zenbuUpload (args, ...)
}

tagError <- function(BED)
  zenbuTag( standardDescription(BED)
          , 'cellomics'
          , libs[bedToSample(BED), "Error"] %>% as.character)

tagRun <- function(BED)
  zenbuTag( standardDescription(BED)
          , 'Run'
          , libs[bedToSample(BED), "Run"] %>% as.character)

tagMeta <- function(BED, TAG)
  zenbuTag( standardDescription(BED)
          , TAG
          , libs[bedToSample(BED), TAG] %>% as.character)
```

Zenbu uploads and tagging
-------------------------

```{r}
samplesToUpload <- rownames(libs) %>%
  sub(RunA, "RunA", .) %>%
  sub(RunB, "RunB", .) %>%
  sub("$", "_R1.bed", .)

sapply(samplesToUpload, safeZupload)
sapply(samplesToUpload, tagError)
sapply(samplesToUpload, tagGroup)
sapply(samplesToUpload, tagMeta, TAG="Run")
sapply(samplesToUpload, tagMeta, TAG="row")
sapply(samplesToUpload, tagMeta, TAG="column")
```
