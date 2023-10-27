import itertools

class Game:
    """
    Class representing a NC(G) game.
    Three games are possible: NC(C_3), NC_00(C_5) and NC_01(C_5)
    """

    def __init__(self, nbPlayers, v0, v1, sym=False):
        assert(nbPlayers == 3 or nbPlayers == 5)
        self.sym = False

        if sym:
            assert(nbPlayers == 5)
            self.sym = sym # Symmetric version -> players receive type 1 as often as type 0

        self.nbPlayers = nbPlayers
        self.questionDistribution = 1/(self.nbPlayers + 1)

        self.v0 = v0
        self.v1 = v1

    def questions(self):
        '''Return a generator over all questions'''
        if self.sym:
            return self.questionsSym()
        else:
            return self.questionsClassic()

    def questionsClassic(self):
        """
        Generator on all question.
        """
        for i in range(self.nbPlayers):
            yield '0' * i + '1' + '0' * (self.nbPlayers - 1 - i)

        yield '1' * self.nbPlayers

    def questionsSym(self):
        '''
        Generate question for the symmetric version
        '''
        for firstOne in range(self.nbPlayers):
            secondOne = (firstOne + 2) % self.nbPlayers #There is two type 1 per question.
            firstOne, secondOne = min(firstOne, secondOne), max(firstOne, secondOne)
            yield '0' * firstOne + '1' + '0' * (secondOne - firstOne - 1) + '1' + '0' * (self.nbPlayers - secondOne - 1)

        yield '1' * self.nbPlayers

    def involvedPlayers(self, question):
        '''
        return the set of involved players for a specific question
        '''

        #If the question is 111 or 11111
        if question.count("1") == self.nbPlayers:
            return list(range(self.nbPlayers))

        playedType1 = question.index('1')
        # Symmetric question, we must choose the correct player with type 1.
        if question.count("1") == 2:
            if question[playedType1 + 2] != "1":
                playedType1 += 3

        involvedPlayers = [(playedType1 - 1) % self.nbPlayers, playedType1, (playedType1 + 1) % self.nbPlayers]

        return involvedPlayers


    def validAnswer(self, answer, question):
        #When every player receives type 1, an answer is correct when the sum of their answer is odd.
        if not '0' in question:
            return sum([int(bit) for bit in answer]) % 2

        #Else, we check if the sum of invovled players is even
        parity = sum([int(answer[idx]) for idx in self.involvedPlayers(question)]) % 2
        return not parity

    def validAnswerIt(self, question):
        for answer in itertools.product(['0', '1'], repeat=self.nbPlayers):
            answer = "".join(answer)
            if self.validAnswer(answer, question):
                yield answer

    def wrongAnswerIt(self, question):
        for answer in itertools.product(['0', '1'], repeat=self.nbPlayers):
            answer = "".join(answer)
            if not self.validAnswer(answer, question):
                yield answer

    def answerPayoutWin(self, answer):
        '''
        Return the mean payout of an answer.
        '''
        nbOfOne = sum((int(bit) for bit in answer))
        return 1/self.nbPlayers * (nbOfOne * self.v1 + (self.nbPlayers - nbOfOne) * self.v0)

    def playerPayoutWin(self, answer, playerId):
        '''
        Return the payout of a specific player.
        '''
        playerAnswer = answer[playerId]
        return self.v1 * (playerAnswer == '1') + self.v0 * (playerAnswer == '0')

    def notPlayerPayoutWin(self, answer, playerId):
        '''
        Return the payout of an answer for a given question if the player i answer: not(adivce)
        '''
        playerAnswer = answer[playerId]
        return self.v0 * (playerAnswer == '1') + self.v1 * (playerAnswer == '0')