
# Genetic Algorithm Based Network Generator
Recently, we have more and more network (mostly undirected binary) generators are proposed to get synthetic network *similar* to target network(s). For instance, ERGM tries to get similar (weighted) graphlets structure, or dk-series tries to capture the distribution of graphlets. Some other more *heuristic* models take specific measures like *degree distribution*, *betweenness* as objectives and optimize the parameters. They get pretty good performance on some networks without too much knowledge. 

At the same time, I am thinking, although these models give pretty reasonable explanation on their parameterization and (sometimes) estimation, are they really indicating such *generating process*? So here I am working on a kinda *null model* for testing

> If we evolve random networks towards some objectives(for now, measures), can we get selected networks similar to the target networks?

Here I used word *selected* because such networks are evolved and selected based on fitness. You can also interpret this process (evolving and selecting) as synthesis of networks with parameter of objectives. That is, we don't parameterize the model, the measured objectives like degree/betweenness distribution are directly our parameters. The synthesis process is done by GA (searching).

### Reason for GA
- Finding networks is close to discrete optimization problems
- GA has a mature framework for multi-objective optimization

### Proposed Steps
- [ ] single target network, single objective (KS of degree distribution)
- [ ] single target network, multi-objective (degree, betweenness and clustering coefficient)
- [ ] run this model on multiple networks (not sure)

### Dataset
In this project, I will use brain functional connectivity networks as dataset. General idea of this processing is introduced in [here][id1].

## Reference
[id1]: https://www.nature.com/articles/nrn2575 (Bullmore, Ed, and Olaf Sporns. "Complex brain networks: graph theoretical analysis of structural and functional systems." Nature Reviews Neuroscience 10.3 (2009): 186.)


```python

```
