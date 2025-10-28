"""
Dynamic array/list utilities for UVAFME model.

This module provides dynamically growing arrays similar to the Fortran
vararray module originally written by Arjen Markus.

The module is straightforward: it defines a suitable data structure,
data can be added to the list and you can retrieve data from it.

Original Fortran module written by Arjen Markus.
Translated to Python for the GoG/UVAFME vegetation model.
"""

from typing import List, Optional, Any, Generic, TypeVar
from copy import deepcopy

T = TypeVar('T')


class ListData:
    """
    Base class for list data items.
    In the original Fortran, this was CHAR_DATA with a string value.
    """
    def __init__(self, value: str = ''):
        self.value = value

    def __str__(self):
        return self.value

    def __repr__(self):
        return f"ListData('{self.value}')"

    def __eq__(self, other):
        if isinstance(other, ListData):
            return self.value == other.value
        elif isinstance(other, str):
            return self.value == other
        return False


# Empty data constant for compatibility
EMPTY_LIST_DATA = ListData('')


class DynamicList(Generic[T]):
    """
    A dynamically growing list that mimics the Fortran vararray behavior.

    This class behaves like a Python list but maintains compatibility
    with the original Fortran interface where items are added and
    retrieved by index.
    """

    def __init__(self, initial_capacity: int = 10):
        """
        Create a new dynamic list.

        Args:
            initial_capacity: Initial capacity of the list
        """
        self._data: List[Optional[T]] = [None] * max(1, initial_capacity)
        self._size = 0
        self._growth_rate = 1.1

    def create(self, capacity: Optional[int] = None):
        """
        Create/recreate the list with specified capacity.

        Args:
            capacity: Initial capacity (optional)
        """
        cap = max(1, capacity) if capacity else 10
        self._data = [None] * cap
        self._size = 0

    def destroy(self):
        """Destroy the list and free memory."""
        self._data = []
        self._size = 0

    def size(self) -> int:
        """Return the number of elements in use."""
        return self._size

    def append(self, data: T):
        """
        Append a value to the list.

        Args:
            data: Data to be appended
        """
        if self._size >= len(self._data):
            self._increase_capacity(self._size + 1)

        self._data[self._size] = data
        self._size += 1

    def at(self, index: int) -> Optional[T]:
        """
        Get the value of the nth element of the list (1-based indexing).

        Args:
            index: Index of the element (1-based, like Fortran)

        Returns:
            The element at the specified index or None if out of bounds
        """
        if index < 1 or index > self._size:
            return None
        return self._data[index - 1]

    def put(self, index: int, data: T):
        """
        Put a value at a specific element of the list (1-based indexing).

        Args:
            index: Index of the element (1-based)
            data: Data to be put in the list
        """
        if index < 1:
            return

        if index > len(self._data):
            self._increase_capacity(index)

        self._size = max(self._size, index)
        self._data[index - 1] = data

    def insert_empty(self, pos: int, number: int = 1):
        """
        Insert one or more empty elements.

        Args:
            pos: Position to insert the empty elements (1-based)
            number: Number of empty elements to insert
        """
        if number < 1 or pos < 1 or pos > self._size:
            return

        if self._size + number >= len(self._data):
            self._increase_capacity(self._size + number)

        # Shift elements to the right
        for i in range(self._size - 1, pos - 2, -1):
            self._data[i + number] = self._data[i]

        # Insert empty elements
        for i in range(number):
            self._data[pos - 1 + i] = None

        self._size += number

    def delete_elements(self, pos: int, number: int = 1):
        """
        Delete one or more elements.

        Args:
            pos: Position to start deletion (1-based)
            number: Number of elements to delete
        """
        if number < 1 or pos < 1 or pos > self._size:
            return

        # Shift elements to the left
        for i in range(pos - 1, self._size - number):
            self._data[i] = self._data[i + number]

        self._size -= number

    def _increase_capacity(self, min_capacity: int):
        """
        Expand the array holding the data.

        Args:
            min_capacity: Minimum required capacity
        """
        new_capacity = max(min_capacity, int(self._growth_rate * len(self._data)))

        if new_capacity > len(self._data):
            new_data = [None] * new_capacity
            # Copy existing data
            for i in range(self._size):
                new_data[i] = self._data[i]
            self._data = new_data

    def to_list(self) -> List[T]:
        """Convert to a Python list containing only the used elements."""
        return [self._data[i] for i in range(self._size) if self._data[i] is not None]

    def __len__(self):
        """Return the number of elements in use."""
        return self._size

    def __getitem__(self, index: int):
        """Get item using 0-based indexing (Python style)."""
        if index < 0 or index >= self._size:
            raise IndexError("Index out of range")
        return self._data[index]

    def __setitem__(self, index: int, value: T):
        """Set item using 0-based indexing (Python style)."""
        if index < 0 or index >= self._size:
            raise IndexError("Index out of range")
        self._data[index] = value

    def __iter__(self):
        """Iterate over the used elements."""
        for i in range(self._size):
            if self._data[i] is not None:
                yield self._data[i]


# Convenience functions that match the Fortran interface
def list_create(capacity: Optional[int] = None) -> DynamicList:
    """Create a new dynamic list."""
    return DynamicList(capacity or 10)


def list_destroy(lst: DynamicList):
    """Destroy a dynamic list."""
    lst.destroy()


def list_size(lst: DynamicList) -> int:
    """Return the number of elements in use."""
    return lst.size()


def list_append(lst: DynamicList, data: Any):
    """Append a value to the list."""
    lst.append(data)


def list_at(lst: DynamicList, index: int) -> Any:
    """Get the value of the nth element (1-based indexing)."""
    return lst.at(index)


def list_put(lst: DynamicList, index: int, data: Any):
    """Put a value at a specific element (1-based indexing)."""
    lst.put(index, data)


def list_insert_empty(lst: DynamicList, pos: int, number: int = 1):
    """Insert empty elements."""
    lst.insert_empty(pos, number)


def list_delete_elements(lst: DynamicList, pos: int, number: int = 1):
    """Delete elements."""
    lst.delete_elements(pos, number)


# Simple list operations for the lists.f90 module
def append_to_species_list(species_list: List[Any], new_species: Any):
    """
    Append a species to a species list.

    This function mimics the append subroutine from lists.f90 that works
    with SpeciesData arrays.
    """
    species_list.append(deepcopy(new_species))


def delete_from_species_list(species_list: List[Any], index: int):
    """
    Delete a species from a species list.

    Args:
        species_list: List of species
        index: Index to delete (1-based, like Fortran)
    """
    if index < 1 or index > len(species_list):
        print(f"List index {index} not in range of list (1-{len(species_list)})")
        return

    species_list.pop(index - 1)  # Convert to 0-based for Python