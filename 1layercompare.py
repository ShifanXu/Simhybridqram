import numpy as np
import qiskit
from qiskit import *
from qiskit import Aer, transpile
from qiskit.tools.visualization import plot_histogram, plot_state_city
import qiskit.quantum_info as qi
from qiskit.providers.aer import PulseSimulator
from qiskit.providers.fake_provider import *
from qiskit_aer.pulse import *
from qiskit.providers import *
from qiskit_experiments.framework import ParallelExperiment
from qiskit_experiments.library import StateTomography
from qiskit_aer import AerSimulator
import qiskit_aer.noise as noise
from qiskit.providers.aer.noise import NoiseModel

'''
ad = QuantumRegister(1, name='ad1')
qram = QuantumRegister(5, name='qram')
bus = QuantumRegister(1, name='bus')
qc = QuantumCircuit(bus,ad,qram)
###address loading
qc.h(ad[0])
qc.swap(ad[0],qram[2])

### retieve data
#qc.barrier()

qc.cnot(qram[2],qram[3])
qc.x(qram[2])
qc.cnot(qram[2],qram[0])
qc.x(qram[2])

#qc.barrier()

qc.swap(qram[0], qram[1])
### qc.swap(qram[4], qram[5])

#qc.barrier()

qc.cnot(qram[1],bus[0])
qc.cnot(qram[4],bus[0])

#qc.barrier()

qc.swap(qram[0], qram[1])

#qc.barrier()

qc.cnot(qram[2],qram[3])
qc.x(qram[2])
qc.cnot(qram[2],qram[0])
qc.x(qram[2])

#qc.barrier()
qc.swap(ad[0],qram[2])
target_state = qi.Statevector.from_instruction(qc)
'''

bus = QuantumRegister(1, name='bus')
ad2 = QuantumRegister(3, name='ad2')
qram = QuantumRegister(16, name='qram')
qc = QuantumCircuit(bus,ad2,qram)
#load address
qc.h(ad2[0])
qc.h(ad2[1])
qc.cnot(ad2[1],ad2[2])
qc.x(ad2[1])
qc.swap(ad2[0],qram[0])
qc.x(qram[0])
qc.cswap(qram[0],ad2[1],qram[1])
qc.cswap(qram[0],ad2[2],qram[2])
qc.x(qram[0])
qc.cswap(qram[0],ad2[1],qram[3])
qc.cswap(qram[0],ad2[2],qram[4])





qc.barrier()

### retieve data

qc.cnot(qram[1],qram[5])
qc.cnot(qram[2],qram[7])
qc.cnot(qram[3],qram[9])
qc.cnot(qram[4],qram[11])
 


qc.barrier()

qc.swap(qram[5],qram[6])
qc.swap(qram[7],qram[8])

qc.barrier()

qc.cnot(qram[6],qram[13])
qc.cnot(qram[8],qram[13])
qc.cnot(qram[10],qram[14])
qc.cnot(qram[12],qram[14])
qc.cnot(qram[13],qram[15])
qc.cnot(qram[14],qram[15])
qc.cnot(qram[15],bus[0])
qc.barrier()

# uncomputing
qc.cnot(qram[13],qram[15])
qc.cnot(qram[14],qram[15])
qc.cnot(qram[6],qram[13])
qc.cnot(qram[8],qram[13])
qc.cnot(qram[10],qram[14])
qc.cnot(qram[12],qram[14])


#Unload address
qc.barrier()
qc.swap(qram[5],qram[6])
qc.swap(qram[7],qram[8])
qc.barrier()

qc.cnot(qram[1],qram[5])
qc.cnot(qram[2],qram[7])
qc.cnot(qram[3],qram[9])
qc.cnot(qram[4],qram[11])
 


qc.barrier()

qc.x(qram[0])
qc.cswap(qram[0],ad2[1],qram[1])
qc.cswap(qram[0],ad2[2],qram[2])
qc.x(qram[0])
qc.cswap(qram[0],ad2[1],qram[3])
qc.cswap(qram[0],ad2[2],qram[4])
qc.swap(ad2[0],qram[0])
qc.x(ad2[1])
qc.cnot(ad2[1],ad2[2])
target_state = qi.Statevector.from_instruction(qc)


# Example error probabilities
from qiskit_aer.noise import (NoiseModel, QuantumError, ReadoutError,
    pauli_error, depolarizing_error, thermal_relaxation_error)
p_gate1 = 0.05

# QuantumError objects

error_gate1 = pauli_error([('Z',p_gate1), ('I', 1 - p_gate1)])
error_gate2 = error_gate1.tensor(error_gate1)

# Add errors to noise model
noise_bit_flip = NoiseModel()
noise_bit_flip.add_all_qubit_quantum_error(error_gate1, ["u1", "u2", "u3"])
noise_bit_flip.add_all_qubit_quantum_error(error_gate2, ["cx"])

print(noise_bit_flip)

F_bell = 0

sim_noise = AerSimulator(noise_model=noise_bit_flip)
qc.save_statevector()
for x in range(1000):
    # Grab results from the job
    result = sim_noise.run(qc).result()
    rho_fit= result.get_statevector(qc)
    F_bell += qi.state_fidelity(rho_fit, target_state)
    #print(qi.state_fidelity(rho_fit, target_state))
F_bell = F_bell / 1001
print('State Fidelity: F = {:.6f}'.format(F_bell))