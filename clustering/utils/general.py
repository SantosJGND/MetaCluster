import os
from random import randint
from typing import Type


def reverse_dict_of_lists(dict_of_lists: dict) -> dict:
    """
    Return dictionary of lists with keys as values and values as keys
    """
    return {value: key for key, values in dict_of_lists.items() for value in values}



class Temp_File:
    """
    Temporary file.
    """

    def __init__(self, temp_dir: str, prefix: str = "temp", suffix: str = ".txt"):
        """
        Initialize.
        """

        self.temp_dir = temp_dir
        self.prefix = prefix
        self.suffix = suffix

        self.path = os.path.join(
            self.temp_dir, f"{self.prefix}_{randint(1000000, 9999999)}{self.suffix}"
        )
        self.file = os.path.basename(self.path)

    def __enter__(self):
        """
        Enter.
        """

        open(self.path, "w").close()
        return self.path

    def __exit__(
        self,
        exc_type: Type[BaseException],
        exc_value: BaseException,
        traceback,
    ):
        """
        Exit.
        """

        if os.path.exists(self.path):
            os.remove(self.path)
