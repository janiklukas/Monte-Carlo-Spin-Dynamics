# Monte Carlo Spin Dynamics

We consider a regular grid of $N^D$ spins described by the Ising model, where $D$ is the number of dimensions. The goal is to calculate the expectation values of the energy and magnetisation using a Markov chain Monte carlo method.

## Ising Model

At each lattice point $\textbf{n}=(n_1,\ldots,n_D)$, the spin is $s(\textbf{n})=\pm 1$. The magnetisation of a spin confuguration $s(\textbf{n})$ is
$$M(s)=\sum_{\textbf{n}}s(\textbf{n}).$$

The energy of a configuration is the sum of internal neighbor interactions and the effect of an external magnetic field $B$:
$$H(s)=-J\sum_{\textbf{n}}\sum_{k=1}^Ds(\textbf{n})s(\textbf{n}+\textbf{e}_k)-\mu B\sum_{\textbf{n}}s(\textbf{n}).$$
We use periodic boundary conditions in evaluating the first term.

The equilibrium probability for a configuration is
$$P(s)=\frac1Z\exp\bigg(-\frac{H(s)}{k_BT}\bigg)\quad\text{with}\quad Z=\sum_s\exp\bigg(-\frac{H(s)}{k_BT}\bigg).$$

## Metropolis Algorithm

To generate a Markov chain of sample spin states, we use the Metropolis algorithm. In each step, it loops over all lattice points $\textbf{n}$ and calculates how the energy would change if that spin was flipped:
$$\Delta H=2J\,s(\textbf{n})\sum_{k=1}^D[s(\textbf{n}+\textbf{e}_k)+s(\textbf{n}-\textbf{e}_k)]+2\mu B\,s(\textbf{n}).$$

The spin $s(\textbf{n})$ is then flipped if $\exp(-\Delta H/k_BT)$ is larger than some random number drawn uniformly from $[0,1]$.

## Covariance Method

The expectation value of an observable $X(s)$ for a given Markov chain of length $N_\text{con}$ is obtained by summing over all thermalized configurations:
$$\bar X=\frac{1}{N_\text{con}-N_\text{th}}\sum_{k=N_\text{th}+1}^{N_\text{con}}X(s_k).$$

The corresponding statistical variance can be estimated using
$$\delta(W)=\frac{1}{(N_\text{con}-N_\text{th})^2}\sum_{m=N_\text{th}}^{N_\text{con}}\sum_{\substack{n=N_\text{th} \\ |m-n|\leq W}}^{N_\text{con}} [X(s_m)-\bar X][X(s_n)-\bar X]$$
where $W$ is a summation window, chosen such that
$$\bigg\vert\frac{\delta(W+1)-\delta(W)}{\delta(W)}\bigg\vert$$
is smaller than some tolerance value to reduce noise in the sum. The standard deviation is then given by $\sqrt{\delta(W)}$.