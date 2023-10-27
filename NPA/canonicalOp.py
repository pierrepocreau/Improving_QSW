from itertools import groupby

def opToPlayer(op, operatorsPlayers):
    """
    Return the player Id with which the operator is associated.
    """
    assert(op != 0) #Id operator associated to every player.
    return (op - 1) // 2

def simplify(monomial, playersOperators, monomialSize=None):
    if monomialSize == None:
        monomialSize = len(monomial)

    canonical = []

    operators = filter(lambda op: op != 0, monomial)  # filter out identity
    operators = sorted(operators, key=lambda op: opToPlayer(op,
                                                            playersOperators))  # Sort without commuting operators of a same player

    for k, g in groupby(operators):
        # GroupBy because the operators are projectors, so op^2 = op
        # if k so that 0 (which correspond to the Id operator) are ignore we put them at the end.
        if k: canonical.append(k)

    # Fill the end with 0 (Id operators)
    while len(canonical) < monomialSize:
        canonical.append(0)

    return canonical

class CanonicalMonomial:

    def __init__(self, monomialList, i, j, playersOperators):
        """
        create the proba/variable <phi| Si^Dagger Sj |phi>
        """
        self.op_i = list(reversed(monomialList[i]))
        self.op_j = monomialList[j]
        self.value = self.op_i + self.op_j

        self.playersOperators = playersOperators #Operators associated to each player
        self.canonical = self.canonicalForm()

    def canonicalForm(self):
        #The moment matrix is symmetric, therefore a monomial and the one where i and j are swapped have the same canonical form.
        return min(simplify(self.value, self.playersOperators), simplify(list(reversed(self.value)), self.playersOperators))

    def __eq__(self, other):
        if not isinstance(other, CanonicalMonomial):
            return False

        #Two operator are equal if they have the same canonic form.
        return other.canonical == self.canonical

    def __hash__(self):
        return tuple(self.canonical).__hash__()