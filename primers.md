Some primer sequences:

````
   PCR                                      TCGTCGGCAGCGTCAGATGTG
    RT                                      TCGTCGGCAGCGTCAGATGTGNNNNNN
   TSO                                      TCGTCGGCAGCGTCAGATGTGTATAAGAGACAGNNNNNNNNTATA(rG)(rG)(rG)
dirFMi AATGATACGGCGACCACCGAGATCTACACiiiiiiiiTCGTCGGCAGCGTCAGATGTGTATAAGAGACAG

N7xx   CAAGCAGAAGACGGCATACGAGATiiiiiiiiGTCTCGTGGGCTCGG
````

The `RT` oligonucleotide primers the reverse-transcription of the mRNA and the 
`TSO` templates the extension of the first-strand cDNA.  A single `PCR` primer
amplifies the cDNAs.  After tagmentation, the library is finally amplified with
indexed primers: `dirFMi` matches the TSO sequence and `N7xx` matches the
tagmentation ends.
