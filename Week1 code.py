# %% Hidden code exerceise 1.1
# all lower text (good practice)
def PatternCount(Text, Pattern):
    """Count the input patern in the text."""
    Text = str(Text)
    Pattern = str(Pattern)
    count = 0
    for i in range(len(Text) - len(Pattern) + 1):
        if Text.find(Pattern, i, i + len(Pattern)) != -1:
            count = count + 1
    return count


# %% Hidden code exercise 1.2
def FrequentWords(text, k):
    """Generate k-frequent words of text."""
    thisdict = {}
    for i in range(len(text) - k + 1):
        kmer = text[i: (i + k)]
        # print(kmer)
        try:
            thisdict[kmer] = thisdict[kmer] + 1
            # print(thisdict.keys())
        except KeyError:
            thisdict.update({kmer: 1})
    maxcount = max(thisdict.values())
    maxlist = []
    for i in thisdict.keys():
        if thisdict[i] == maxcount:
            maxlist.append(i)
    return maxlist


# %% Hideen code exercise 1.3
def ReverseComplement(Pattern):
    """Find the reverse complement of a DNA string."""
    thisdict = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    complement = []
    for i in Pattern:
        complement.append(thisdict[i])
    return "".join(complement)[::-1]


# %% Hidden code exercise 1.4
def PatternMatching(Pattern, Genome):
    """Find all occurrences of a pattern in a string."""
    occurrences = []
    for i in range(len(Genome)):
        if Pattern == Genome[i: (i + len(Pattern))]:
            occurrences.append(str(i))
    return ' '.join(occurrences)


# %% Hidden code exercise 1.5
def ClumpFinding(genome, k, L, t):
    """Find patterns forming clumps in a string."""
    clumps = set()
    for i in range(len(genome)):
        text = genome[i: (i + L)]
        thisdict = {}
        for j in range(L):
            kmers = text[j: (j + k)]
            try:
                thisdict[kmers] = thisdict[kmers] + 1
            except KeyError:
                thisdict.update({kmers: 1})
        for m in thisdict.keys():
            if thisdict[m] == t:
                clumps.add(m)
    return clumps


# %% Hidden code exercise 1.6
_thisdict_ = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}


def patterntonumber(text, k):
    """Assign strings to numbers."""
    # print(text)
    count = k - 1
    number = 0
    for i in text:
        if i == 'A':
            number = number + 4 ** count * 0
        elif i == 'C':
            number = number + 4 ** count * 1
        elif i == 'G':
            number = number + 4 ** count * 2
        elif i == 'T':
            number = number + 4 ** count * 3
        count -= 1
    # print(number)
    return number


def ComputingFrequencies(text, k):
    """Generate a frequency array."""
    frequency = []
    for _i in range(4 ** k):
        frequency.append(0)
    for i in range(len(text) - k + 1):
        frequency[patterntonumber(text[i: (i + k)], k)] += 1
    return frequency


# %% Hidden code excercise 1.8
def NumberToPattern(index, k):
    """Transform numbers to nucleotide codons."""
    thislist = []
    thisdict = {0: 'A', 1: 'C', 2: 'G', 3: 'T'}
    if index <= 4 ** k:
        while k > 0:
            thislist.append(thisdict[int(index / 4 ** (k - 1))])
            index = index % (4 ** (k - 1))
            k -= 1
        return "".join(thislist)
    else:
        return print('Please enter a valid numbers. Index must be small than 4^k')
