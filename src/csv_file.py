"""
CSV file utilities for writing structured data.

This module provides functionality for writing CSV files in a format
compatible with the original Fortran UVAFME model.

Original Fortran module written by Arjen Markus.
Translated to Python for the GoG/UVAFME vegetation model.
"""

import csv
import io
from typing import Union, List, Any, Optional, TextIO
import numpy as np


class CSVWriter:
    """
    A CSV writer class that mimics the behavior of the Fortran csv_file module.

    This class maintains compatibility with the original Fortran interface
    where file units are used and comma-separated values are written with
    specific formatting rules.
    """

    def __init__(self, file_handle: TextIO):
        """
        Initialize CSV writer with a file handle.

        Args:
            file_handle: Open file handle for writing
        """
        self.file_handle = file_handle
        self.current_record = []

    def csv_next_record(self):
        """
        Go to the next record (convenience routine).

        The current record is closed, the next write will be to the new record.
        """
        if self.current_record:
            self.file_handle.write(','.join(self.current_record) + '\n')
            self.current_record = []
        else:
            self.file_handle.write('\n')

    def csv_write(self, value: Any, advance: bool = True):
        """
        Generic CSV write function that dispatches to specific type handlers.

        Args:
            value: Value to write (int, float, str, or array)
            advance: Whether to advance to next record after writing
        """
        if isinstance(value, (list, np.ndarray)):
            if hasattr(value, 'ndim') and value.ndim == 2:
                self._csv_write_2d(value)
            else:
                self._csv_write_1d(value, advance)
        elif isinstance(value, int):
            self._csv_write_integer(value, advance)
        elif isinstance(value, float):
            self._csv_write_real(value, advance)
        elif isinstance(value, str):
            self._csv_write_char(value, advance)
        else:
            # Try to convert to string as fallback
            self._csv_write_char(str(value), advance)

    def _csv_write_integer(self, value: int, advance: bool):
        """Write a single integer to the CSV file."""
        formatted_value = f"{value:d}"
        self.current_record.append(formatted_value)

        if advance:
            self.csv_next_record()

    def _csv_write_real(self, value: float, advance: bool):
        """Write a single real (float) to the CSV file."""
        # Use scientific notation similar to Fortran G format
        if abs(value) < 1e-6 or abs(value) >= 1e6:
            formatted_value = f"{value:.6e}"
        else:
            formatted_value = f"{value:.6f}"
            # Only strip trailing zeros for non-scientific notation
            formatted_value = formatted_value.rstrip('0').rstrip('.')

        self.current_record.append(formatted_value)

        if advance:
            self.csv_next_record()

    def _csv_write_char(self, value: str, advance: bool):
        """Write a single character string to the CSV file."""
        # Handle embedded quotes by doubling them
        escaped_value = value.replace('"', '""')

        # Quote the string if it contains commas, quotes, or newlines
        if ',' in escaped_value or '"' in value or '\n' in escaped_value:
            formatted_value = f'"{escaped_value}"'
        else:
            formatted_value = escaped_value.strip()

        self.current_record.append(formatted_value)

        if advance:
            self.csv_next_record()

    def _csv_write_1d(self, array: Union[List, np.ndarray], advance: bool = True):
        """Write a one-dimensional array to the CSV file."""
        for i, item in enumerate(array):
            is_last = (i == len(array) - 1)
            self.csv_write(item, advance=(is_last and advance))

    def _csv_write_2d(self, array: Union[List[List], np.ndarray]):
        """Write a two-dimensional array to the CSV file."""
        if hasattr(array, 'shape'):
            # NumPy array
            for row in array:
                self._csv_write_1d(row, advance=True)
        else:
            # List of lists
            for row in array:
                self._csv_write_1d(row, advance=True)


# Global dictionary to store file handles and their associated CSV writers
_csv_writers = {}


def csv_write(file_unit: Union[int, TextIO], value: Any, advance: bool = True):
    """
    Write data to a CSV file using a file unit or file handle.

    This function maintains compatibility with the Fortran interface where
    file units (integers) are used to identify open files.

    Args:
        file_unit: File unit number (int) or file handle (TextIO)
        value: Value to write
        advance: Whether to advance to next record
    """
    # Get or create CSV writer for this file unit
    if file_unit not in _csv_writers:
        if isinstance(file_unit, int):
            # For integer file units, we assume the file handle is managed elsewhere
            # This is a limitation of the Python translation - in practice,
            # file handles should be passed directly
            raise ValueError(f"File unit {file_unit} not found. Pass file handle directly.")
        else:
            _csv_writers[file_unit] = CSVWriter(file_unit)

    writer = _csv_writers[file_unit]
    writer.csv_write(value, advance)


def csv_next_record(file_unit: Union[int, TextIO]):
    """
    Advance to the next record in the CSV file.

    Args:
        file_unit: File unit number or file handle
    """
    if file_unit in _csv_writers:
        _csv_writers[file_unit].csv_next_record()


def register_csv_file(file_unit: int, file_handle: TextIO):
    """
    Register a file handle with a file unit number for CSV writing.

    Args:
        file_unit: File unit number
        file_handle: Open file handle
    """
    _csv_writers[file_unit] = CSVWriter(file_handle)


def close_csv_file(file_unit: Union[int, TextIO]):
    """
    Close and cleanup CSV writer for a file unit.

    Args:
        file_unit: File unit number or file handle
    """
    if file_unit in _csv_writers:
        # Flush any remaining data
        _csv_writers[file_unit].csv_next_record()
        del _csv_writers[file_unit]


# Convenience functions for specific data types (maintaining Fortran compatibility)
def csv_write_integer(file_unit: Union[int, TextIO], value: int, advance: bool = True):
    """Write a single integer to CSV file."""
    csv_write(file_unit, value, advance)


def csv_write_real(file_unit: Union[int, TextIO], value: float, advance: bool = True):
    """Write a single float to CSV file."""
    csv_write(file_unit, value, advance)


def csv_write_char(file_unit: Union[int, TextIO], value: str, advance: bool = True):
    """Write a single string to CSV file."""
    csv_write(file_unit, value, advance)


def csv_write_integer_1d(file_unit: Union[int, TextIO], array: List[int], advance: bool = True):
    """Write a 1D array of integers to CSV file."""
    csv_write(file_unit, array, advance)


def csv_write_real_1d(file_unit: Union[int, TextIO], array: List[float], advance: bool = True):
    """Write a 1D array of floats to CSV file."""
    csv_write(file_unit, array, advance)


def csv_write_char_1d(file_unit: Union[int, TextIO], array: List[str], advance: bool = True):
    """Write a 1D array of strings to CSV file."""
    csv_write(file_unit, array, advance)


def csv_write_integer_2d(file_unit: Union[int, TextIO], array: List[List[int]]):
    """Write a 2D array of integers to CSV file."""
    csv_write(file_unit, array)


def csv_write_real_2d(file_unit: Union[int, TextIO], array: List[List[float]]):
    """Write a 2D array of floats to CSV file."""
    csv_write(file_unit, array)


def csv_write_char_2d(file_unit: Union[int, TextIO], array: List[List[str]]):
    """Write a 2D array of strings to CSV file."""
    csv_write(file_unit, array)