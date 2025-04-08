# Lyapunov Functions in Model Predictive Control (MPC)

Dear Majid,  
you raise an interesting question.

Let's consider a discrete-time Optimal Control Problem (OCP) which implicitly defines an MPC feedback:

$$
V_N(\hat x_i) := \min \sum_{k=0}^{N-1} l(x_k, u_k) + V_f(x_N)
$$

**Subject to:**

- \( x_{k+1} = f(x_k, u_k), \quad x_0 = \hat x_i \)
- \( x_k \in \mathbb{X}, \quad u_k \in \mathbb{U} \quad \forall k \in \{0, ..., N\} \)
- \( x_N \in \mathbb{X}_f \)

Here:

- \( l(x_k, u_k) \) is the cost function or **stage cost**
- \( V_f(x_N) \) is the **terminal cost**
- \( V_N(\hat x_i) \) is referred to as the **optimal value function** of the OCP

---

For most MPC designs (e.g., set stabilization, tracking), one proves stability by choosing \( V_N(\hat x_i) \) as a **Lyapunov function candidate**. See, for example, Mayne et al. (2000) or Chen & Allgöwer (1998).

In such proofs, it's typically required that the terminal cost \( V_f(x_N) \) satisfies a **local Lyapunov-like decrease property** on the terminal set \( \mathbb{X}_f \) to ensure **asymptotic stability**. That is, it is required that:

$$
\forall x \in \mathbb{X}_f, \exists u \in \mathbb{U} \text{ such that } V_f(f(x,u)) - V_f(x) \leq -l(x, u)
$$

---

To come back to your question: the typical Lyapunov function candidate in MPC is the **optimal value function** \( V_N \). However, to show that \( V_N \) is a Lyapunov function of the closed loop, many works require that the **terminal cost** satisfies a **Lyapunov-like decrease property** on the terminal set.

Other works directly consider \( V_f \) to be a **control Lyapunov function**, see Jadbabaie et al. (2001). Also, note that in **economic MPC**, one often considers **shifted value functions** as Lyapunov function candidates.

I hope this helps to understand what kind of role Lyapunov functions play in the stability of MPC.

---

## References

- Mayne, David Q., et al. *"Constrained model predictive control: Stability and optimality."* Automatica 36.6 (2000): 789-814.
- Chen, H., and F. Allgöwer. *"A Quasi-Infinite Horizon Nonlinear Model Predictive Control Scheme with Guaranteed Stability."* Automatica 34.10 (1998): 1205-1217.
- Jadbabaie, Ali, Jie Yu, and John Hauser. *"Unconstrained receding-horizon control of nonlinear systems."* IEEE Transactions on Automatic Control 46.5 (2001): 776-783.
