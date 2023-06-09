# Notes for Talk

## Jun
### Slide 1 (Title): 
Today, we want to present on our progress for our Qiskit Advocate Mentorship project.

### Slide 2 (Paper #1):
So for the past 3 months, we have been studying in depth two papers:

1. Constructive Quantum Shannon Decomposition from Cartan Involutions by Byron Drury and Peter Love

### Slide 3 (What is KAK?):
This paper describes a general algorithm for decomposing circuit using a theorem called Cartan KAK decomposition. However for our project, we focus on the case n=2, that is to say two qubit decomposition. This case serves as a base case in more general routine like Quantum Shannon Decomposition described in the paper.

In addition to trying to understand their results, we also produce a set of notes to guide future readers through the technicalities. We also implement the KAK decomposition in python and created a detailed tutorial on our process.

## Minh

### Slide 4 (Paper #2):
The second paper we studied is 
2. A geometric theory of non-local two-qubit operations by Zhang et al.

This is in some sense an extension of the first paper. For quantum with one qubit, we have superposition. But with two qubit, we have entanglement. So this paper outlines a geometric method for study two-qubit entanglement. By using tools from Lie theory, the author showed that we can house all the two qubit gates under a Weyl chamber in the shape of a tetrahedron that lives in 3-d space. More specifically, if we’re given a two qubit gate, we can find the unique so called canonical coordinates for the gate in the Well chamber.

### Slide 5 (Weyl Chamber):
Now then, we can meaningfully ask question like: What proportion two qubit gates can produce the bell state? And so when we look at the structure of this set of bell-state producing gates they form another geometric object inside the Weyl chamber, and now all we need to do calculate the volume of this object and divide by the volume of the tetrahedron. The answer actually turns out to be 1/2. That is 1/2 of two qubit gates produces bell state.

So for this paper, we also wrote up some notes and implement the codes to calculate the canonical coordinates.

### Slide 6 (Future Directions):
So going forward, we are working on a medium article to give a high level overview of these two papers. The article will give a little more motivation and more visualization to help the readers understand these somewhat technical results.

### Slide 7 (Links):
Here are the links to our notes and code implementation
