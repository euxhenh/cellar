from abc import ABC, abstractmethod


class Unit(ABC):
    def __init__(self):
        pass

    @abstractmethod
    def get(self):
        pass
