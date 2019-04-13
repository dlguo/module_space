# Research repository

My temporary unorganized research project and experiment. I will remove the content from here and create separate project/repository after the project is completed. Current progress:

- [x] Entropy measure for network variability
- [x] Variational auto-encoder for vectorized networks
- [ ] TSP reordered network CNN encoder
- [ ] Differentiable variability loss function

## Dataset

HCP data has 7 different tasks and 2 separate resting-state scan sessions. The protocol is as follows:

## Data description
All data used here is from HCP 900 subjects. There are 100 randomly picked subjects in this dataset and 2 resting-state scan sessions (4 files) and 7 task sessions (14 files) are done w/ them. The protocol are listed as below:

Condition | Runs | Frames per run | Run Duration (min:sec) 
--- | --- | --- | ---
REST (Resting-state) | 4 | 1200 | 14:33 
Working Memory | 2 | 405 | 5:01
Gambling | 2 | 253 | 3:12
Motor | 2 | 284 | 3:34
Language | 2 | 316 | 3:57
Social Cognition | 2 | 274 | 3:27
Relational Processing | 2 | 232 | 2:56
Emotion Processing | 2 | 176 | 2:16

To exclude the impact of session length in this experiment(samples v.s. correlation), I randomly select the time series of 150-step (a bit higher than the shortest session, emotion) and then do the correlation. The correlation/adjacency matrices used later are all processed as this unless stated otherwise.