# Tasks ToDo
- [x] download ToxProt 03.2025 (aka `Release 2025_01`)
- differentiate between marine and terrestrial (inc. freshwater)
  - [x] get family level taxonomy
  - [x] Ivan manually classifies them

## figures
- [ ] table comparison (2017 to 2025)
  - [x] number protein entries
  - [x] unique species
  - [x] number taxa - aka higher rank (order / family / genus)
  - [x] protein families (regex to improve)
  - [x] fragments
  - [x] PTM annotation
  - [x] toxic dose
- comparison (2017 - 2025)
  - [x] protein family of top 10 + other + NaN  - stackedbar plot
  - [x] protein length distribution (25 steps and final is 300+) - overlayed histogram
  - [x] taxa distribution of top 5 + other - stacked bar plot
  - [x] taxa newcomer (optional)
  - [x] terrestrial vs marine protein family 17 vs 25 - grouped bar plot
- protspace plot of final dataset
  - [ ] new regex
  - [ ] remove fragments + SP

## other
- color schema - tab10

## Task distribution
- ST
  - [x] protein families clustering
  - [x] createtable comparison - matplotlib?
  - [x] create plots (protein length distribution + protein family)
  - [ ] apply new regex/method to compress/summarize taxa
  - [ ] Count by spiecies instead of scientific name
- TS
  - [x] Parse kim excel file to CSV into processed
  - [x] Parse original 2017.11 swissprot data
  - [x] taxa extraction (order, family genus) -> taxaiq?
  - [x] distinguish between marine and terrestrial taxa family
  - [x] create plots (tax dist + taxa new + terra vs marine)
  - [ ] ProtSpace plot
  - [ ] Add column "species" instead of "scientific name"



## ToDo
- [ ] fix x-axis on histogram 1-25, 26-50, ... 300+
- [ ] fix protein family plot
  - [ ] Huwenotoxin - neurotoxin 10
  - [ ] sort by abundancy
  - [ ] sankey diagram would be great
- [ ] double check the protein family count
- [ ] Marine/Terrestrial protein families improve plot
- [ ] stacked bar (top 5 taxa) -> dots in background
- [ ] write methods
- [ ] supplements/appendix the download version
- [ ] Gigascience