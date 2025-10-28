"""
Utility functions for UVAFME vegetation model.
Translated from Utilities.f90
"""

import math
import numpy as np
from typing import List, Tuple, Optional


def kron(x: float) -> float:
    """Kronecker delta function.
    Returns 1.0 if x >= 0.0, else 0.0
    """
    return 1.0 if x >= 0.0 else 0.0


def roundtoN(x: float, N: int) -> float:
    """Round x to N decimal places."""
    if N <= 0:
        return round(x)
    else:
        factor = 10.0 ** N
        return round(x * factor) / factor


def stddev(array: List[float], missing: float = -999.0) -> Tuple[float, float]:
    """Compute arithmetic mean and standard deviation of array.
    
    Args:
        array: List of values
        missing: Value to treat as missing data
        
    Returns:
        Tuple of (arithmetic_mean, standard_deviation)
    """
    # Filter out missing values
    valid_values = [x for x in array if x != missing]
    
    if len(valid_values) == 0:
        return 0.0, 0.0
    
    # Calculate mean
    arith_mean = sum(valid_values) / len(valid_values)
    
    if len(valid_values) == 1:
        return arith_mean, 0.0
    
    # Calculate standard deviation
    variance = sum((x - arith_mean) ** 2 for x in valid_values) / (len(valid_values) - 1)
    # Protect against negative variance due to floating point errors
    variance = max(0.0, variance)
    std_deviation = math.sqrt(variance)
    
    return arith_mean, std_deviation


def sort_strings(array: List[str]) -> List[str]:
    """Shell sort for character arrays.
    
    Returns a sorted copy of the input array.
    """
    arr = array.copy()
    length = len(arr)
    
    # Shell sort implementation
    gap = length // 2
    
    while gap > 0:
        for i in range(gap, length):
            temp = arr[i]
            j = i
            
            while j >= gap and arr[j - gap] > temp:
                arr[j] = arr[j - gap]
                j -= gap
            
            arr[j] = temp
        
        gap //= 2
    
    return arr


def safe_divide(numerator: float, denominator: float, default: float = 0.0) -> float:
    """Safely divide two numbers, returning default if denominator is zero."""
    if abs(denominator) < 1e-10:
        return default
    return numerator / denominator


def clamp(value: float, min_val: float, max_val: float) -> float:
    """Clamp value between min_val and max_val."""
    return max(min_val, min(max_val, value))


def interpolate_linear(x: float, x1: float, y1: float, x2: float, y2: float) -> float:
    """Linear interpolation between two points."""
    if abs(x2 - x1) < 1e-10:
        return y1
    
    t = (x - x1) / (x2 - x1)
    return y1 + t * (y2 - y1)


def moving_average(data: List[float], window_size: int) -> List[float]:
    """Calculate moving average with given window size."""
    if window_size <= 0 or window_size > len(data):
        return data.copy()
    
    result = []
    for i in range(len(data)):
        start = max(0, i - window_size // 2)
        end = min(len(data), i + window_size // 2 + 1)
        avg = sum(data[start:end]) / (end - start)
        result.append(avg)
    
    return result


def find_index(array: List[str], target: str) -> int:
    """Find index of target string in array. Returns -1 if not found."""
    try:
        return array.index(target)
    except ValueError:
        return -1


def remove_duplicates(array: List[str]) -> List[str]:
    """Remove duplicate strings while preserving order."""
    seen = set()
    result = []
    for item in array:
        if item not in seen:
            seen.add(item)
            result.append(item)
    return result