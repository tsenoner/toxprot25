# Tasks ToDo
- [x] download ToxProt 03.2025 (aka `Release 2025_01`)
- differentiate between marine and terrestrial (inc. freshwater)
  - [ ] get family level taxonomy
  - [ ] Ivan manually classifies them

## figures
- [ ] table comparison (2017 to 2025)
  - [ ] number protein entries
  - [ ] unique species
  - [ ] number taxa - aka higher rank (order / family / genus)
  - [ ] protein families (regex to improve)
  - [ ] fragments
  - [ ] PTM annotation
  - [ ] toxic dose
- comparison (2017 - 2025)
  - [ ] protein family of top 10 + other + NaN  - stackedbar plot
  - [ ] protein length distribution (25 steps and final is 300+) - overlayed histogram
  - [ ] taxa distribution of top 5 + other - stacked bar plot
  - [ ] taxa newcomer (optional)
  - [ ] terrestrial vs marine protein family 17 vs 25 - grouped bar plot
- protspace plot of final dataset
  - [ ] new regex
  - [ ] remove fragments + SP

## other
- color schema - tab10

## Task distribution
- ST
  - [ ] protein families clustering
  - [ ] createtable comparison - matplotlib?
  - [x] create plots (protein length distribution + protein family)
- TS
  - [ ] Parse kim excel file to CSV into processed
  - [ ] taxa extraction (order, family genus) -> taxaiq?
  - [ ] distinguish between marine and terrestrial taxa family
  - [ ] create plots (tax dist + taxa new + terra vs marin)
  - [ ] ProtSpace plot