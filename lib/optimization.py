from abc import ABCMeta, abstractmethod
import time
import datetime

class Optimizer:
    __metaclass__ = ABCMeta
    @abstractmethod
    def __init__(self, **kwargs):
        pass

    @abstractmethod
    def optimize_embedding(self, embedding_model, loss_function, pair_genetator):
        pass

class SGD(Optimizer):
    def __init__(self, learning_rate=0.1, n_epoch=100, verbose=True):
        self.learning_rate = learning_rate
        self.n_epoch = n_epoch
        self.verbose = verbose

    def optimize_embedding(self, embedding_model, loss_function, pair_generator):
        for epoch in range(self.n_epoch):
            pair_generator.shuffle()
            if self.verbose:
                print "Epoch {} / {} ...".format(epoch+1, self.n_epoch)
            start = time.time()
            for batch in pair_generator():
                x = embedding_model.get_state_vector()
                grad = loss_function.loss_gradient(batch, embedding_model.get_distance_info(batch))
                x -= grad * self.learning_rate
                embedding_model.set_state_vector(x)
            finish = time.time()
            if self.verbose:
                print "Elapsed time: {}s".format(datetime.timedelta(seconds=finish-start))
