import numpy as np

def read_fasta_file(filename):
    """
    Reads the given FASTA file f and returns a dictionary of sequences.

    Lines starting with ';' in the FASTA file are ignored.
    """
    sequences_lines = {}
    current_sequence_lines = None
    with open(filename) as fp:
        for line in fp:
            line = line.strip()
            if line.startswith(';') or not line:
                continue
            if line.startswith('>'):
                sequence_name = line.lstrip('>')
                current_sequence_lines = []
                sequences_lines[sequence_name] = current_sequence_lines
            else:
                if current_sequence_lines is not None:
                    current_sequence_lines.append(line)
    sequences = {}
    for name, lines in sequences_lines.items():
        sequences[name] = ''.join(lines)
    return sequences


class hmm:
    def __init__(self, init_probs, trans_probs, emission_probs):
        self.init_probs = init_probs
        self.trans_probs = trans_probs
        self.emission_probs = emission_probs

        
def translate_indices_to_observations(indices):
    mapping = ['a', 'c', 'g', 't']
    return ''.join(mapping[idx] for idx in indices)


def translate_path_to_indices(path):
    return list(map(lambda x: int(x), path))


def translate_indices_to_path(indices):
    return ''.join([str(i) for i in indices])


def translate_observations_to_indices(obs):
    mapping = {'a': 0, 'c': 1, 'g': 2, 't': 3}
    return [mapping[symbol.lower()] for symbol in obs]


def translate_sequence_to_states(sequence):
    N = len(sequence)
    states = np.array([])
    i = 0

    while i < N:
        nextS, lenA = checkStart(sequence[i: i + 3])
        states = np.append(states, nextS, axis = 0)
        i += lenA
        if states[-1] == 3 or states[-1] == 6 or states[-1] == 9:
            while states[-1] != 15 and states[-1] != 18 and states[-1] != 21:
                states = np.append(states, checkEndF(sequence[i : i + 3]), axis = 0)
                i += 3

        if states[-1] == 24 or states[-1] == 27 or states[-1] == 30:        
            while states[-1] != 36 and states[-1] != 39 and states[-1] != 42:
                states = np.append(states, checkEndR(sequence[i : i + 3]), axis = 0)
                i += 3

    return states


### TEST FOR HMM 7 ###
def createZ7(annotation):
    N = len(annotation)
    i = 0
    Z = np.zeros(N)

    while i < N:
        if i == 0:
            Z[i] = 3
            i += 1
        while annotation[i: i + 3] == 'CCC':
            Z[i: i + 3] = np.array([4, 5, 6])
            i += 3
        while annotation[i: i + 3] == 'RRR':
            Z[i: i + 3] = np.array([2, 1, 0])
            i += 3
        Z[i] = 3
        i += 1
        
    return Z


def createA7(Z):
    A = np.zeros((7, 7))
    for i in range(Z.shape[0] - 1):
        a, b = int(Z[i]), int(Z[i + 1])
        A[a, b] += 1

    for i in range(7):
        A[i] /= np.sum(A[i])

    return A


def createPi7():
    Pi = np.zeros(7)
    Pi[3] = 1

    return Pi


def createPhi7(Z, sequence):
    Phi = np.zeros((7, 4))
    for i in range(Z.shape[0]):
        state = int(Z[i])
        emission = int(sequence[i])
        Phi[state, emission] += 1

    for i in range(7):
        Phi[i] /= np.sum(Phi[i])

    return Phi

### END TEST FOR HMM 7 ###


def viterbi(A, Phi, Pi, sequence):
    N = len(sequence) # Number of steps in the markov chain
    K = 7 # Number of hidden states
    Omega = np.zeros((K, N))
    OmegaBack = np.zeros((K, N))

    # First column
    for i in range(K):
        Omega[i, 0] = Pi[i] * Phi[i, sequence[0]]

    # Probably need log to make this work
    for i in range(1, N): # Loop over the sequence
        for j in range(K): # Loop over the hidden states
            preMax = 0
            argpreMax = 0
            for k in range(K): # Loop over previous column (maybe use argmax)
                newPreMax = Omega[k, i - 1] * A[k, j] * Phi[j, sequence[i]]
                if newPreMax > preMax:
                    preMax = newPreMax
                    argpreMax = k
            Omega[j, i] = preMax
            OmegaBack[j, i] = argpreMax

    # Now find the way back
    Z = np.zeros(N)
    y = np.argmax(Omega[:, N - 1]) # Last column 
    Z[-1] = y

    for i in range(N, 1):
        y = OmegaBack[y, i]
        Z[N - 1] = y

    return Z


def translate_sequence_to_states2(sequence, annotation):
    N = len(sequence)
    states = np.array([])
    i = 0
    while i < N:
        if annotation[i: i + 3] == 'CCC' and isStartF(sequence[i: i + 3]):
            states = np.append(states, checkStart(sequence[i: i + 3])[0])
            while annotation[i: i + 3] == 'CCC' or not isStopF(sequence[i : i + 3]):
                states = np.append(states, [10, 11, 12], axis = 0)
                i += 3
            states = np.append(states, checkEndF(sequence[i : i + 3]), axis = 0)
            i += 1
            
        elif annotation[i:i + 3] == 'RRR' and isStartR(sequence[i:i+3]):
            if states[-1] == 24 or states[-1] == 27 or states[-1] == 30:        
                while annotation[i + 1] == 'R' or not isStopR(sequence[i : i + 3]):
                    states = np.append(states, [31, 32, 33], axis = 0)
                    i += 3
                states = np.append(states, checkEndR(sequence[i : i + 3]), axis = 0)
                i += 1
        else:
            states = np.append(states, [0])
            i += 1

        
def isStartF(s):
    if s == 'ATG' or s == 'GTG' or s == 'TTG':
        return True
    else:
        return False

    
def isStartR(s):
    if s == 'TTA' or s == 'CTA' or s == 'TCA':
        return True
    else:
        return False

def isStopF(s):
    if s == 'TAG' or s == 'TGA' or s == 'TAA':
        return True
    else:
        return False

def isStopR(s):
    if s == 'CAT' or s == 'CAC' or s == 'CAA':
        return True
    else:
        return False


def checkStart(string):
    if string == 'ATG':
        return np.array([1, 2, 3]), 3

    if string == 'GTG':
        return np.array([4, 5, 6]), 3

    if string == 'TTG':
        return np.array([7, 8, 9]), 3

    if string == 'TTA':
        return np.array([22, 23, 24]), 3

    if string == 'CTA':
        return np.array([25, 26, 27]), 3

    if string == 'TCA':
        return np.array([28, 29, 30]), 3 

    return np.array([0]), 1


def checkEndF(string):
    if string == 'TAG':
        return np.array([13, 14, 15])

    if string == 'TGA':
        return np.array([16, 17, 18])

    if string == 'TAA':
        return np.array([19, 20, 21])

    return np.array([10, 11, 12])


def checkEndR(string):
    if string == 'CAT':
        return np.array([34, 35, 36])

    if string == 'CAC':
        return np.array([37, 38, 39])

    if string == 'CAA':
        return np.array([40, 41, 42])

    return np.array([31, 32, 33])


def calculateA(states):
    A = np.zeros((42, 42))
    for i in range(states.shape[0]-1):
        a, b = states[i], states[i + 1]
        A[a, b] += 1

    for i in range(42):
        A[i] /= np.sum(A[i])

    return A


def calculatePi():
    pi = np.zeros(42)
    pi[0] = 4
    pi[7] = 1
        


genomes = {}
for i in range(1, 11):
    sequence = read_fasta_file('genome' + str(i) + '.fa')
    genomes['genome' + str(i)] = sequence['genome' + str(i)]
annotation = read_fasta_file('true-ann1.fa')

# Test for hmm7
Z = createZ7(annotation['true-ann1'])
A = createA7(Z)
sequence = translate_observations_to_indices(genomes['genome1'])
Phi = createPhi7(Z, sequence)
Pi = createPi7()
print('Transition probabilities are', A)
print('Emission probabilities are', Phi)
Zml = viterbi(A, Phi, Pi, sequence)
print(Zml)


#states = translate_sequence_to_states(genomes['genome1'])
#np.save('genome1.npy', states)
