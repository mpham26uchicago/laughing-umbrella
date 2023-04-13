1. **Optimal synthesis into fixed XX interactions**: This paper discusses the scenario where you have access to not only CNOT but also 1/2 CNOT, 1/3 CNOT. With this set of 2q gates, you can achieve almost optimal (optimal = having all fraction of the CNOT gate) circuit length. Notice that for our case, we have all fraction of the Molmer-Sorensen gate...
- Github: https://qiskit.org/documentation/_modules/qiskit/quantum_info/synthesis/xx_decompose/decomposer.html
- Talk given by Eric P: https://chromotopy.org/latex/talks/xx-synthesis.pdf, video: https://www.youtube.com/watch?v=8zB2LjUR2FU

2. **Quantum Instruction Set Design for Performance**: This paper discusses the sqrt-iSWAP gate (they call it SQiSW). Amazingly, they show that with just 2-applications of SQiSW gate, you can synthesize 69% of all the 2-qubit unitaries. This paper might be slightly less relevant since we have much more than just 1 gate, but a whole family of gates, but the techniques used should be super helpful.
- Github: https://quantumai.google/reference/python/cirq/two_qubit_matrix_to_sqrt_iswap_operations

3. **Designing calibration and expressivity-efficient instruction sets for quantum computing**: The above 2 papers are analytical, thus are optimal and super efficient. In the case where we cannot find analytical solutions, we will have to be similar to this paper and perform numerical exhaustive search.

4. **Fixed-Depth Two-Qubit Circuits and the Monodromy Polytope**: This paper provides the theory of paper 1 & 2. You will probably need to understand at least part of it.
- ppt: https://chromotopy.org/latex/talks/iwqc19.pdf
- code: https://github.com/Qiskit-Extensions/monodromy
