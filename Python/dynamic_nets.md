
# Dynamic brain networks
This project currently is only a experimental project for finding interesting patterns of dynamic brain functional networks. The data used here is **dynamic FC series** therefore all analysis is based on the dynamics of brain network series. The network series is generated from the Pearson correlation between 360 ROIs (Glasser's parcellation)

Here are some targets for this research:
- [ ] Basic FC series dynamics
- [ ] Modular dynamics
- [ ] Entropy analysis

## Data description
All data used here is from HCP 900 subjects. There are 100 randomly picked subjects in this dataset and 2 resting-state scan sessions(4 fMRI time courses) and 7 task sessions(14 fMRI) are done w/ them. The protocol are listed as below:

Condition | Runs | Frames per run | Run Duration (min:sec) 
--- | --- | ---
REST (Resting-state) | 4 | 1200 | 14:33 
Working Memory | 2 | 405 | 5:01
Gambling | 2 | 253 | 3:12
Motor | 2 | 284 | 3:34
Language | 2 | 316 | 3:57
Social Cognition | 2 | 274 | 3:27
Relational Processing | 2 | 232 | 2:56
Emotion Processing | 2 | 176 | 2:16

## Prepare the data
Because network series data is very huge, we can only do experiments on small section of subjects. Here we picked 4 subjects and try to learn from these data.


```R
# 
subj_list <- c("123925", "125525", "122620", "245333")
rsn_names <- c('visual', 'motor', 'da', 'va', 'limbic', 'fp', 'DMN', 'sub')
sess_list <- list.files('/data/results_SIFT2/120111/fMRI/')[11:12]
```

    [1] 12

