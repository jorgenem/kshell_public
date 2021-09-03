import pytest
from kshell_ui import split_jpn

def test_split_jpn_value_error():
    """
    Test that ValueError is raised if jpn does not contain parity
    information ('+' or '-').
    """
    with pytest.raises(ValueError):
        split_jpn(jpn="", valence_p_n=[0, 0])

def test_split_jpn():
    """
    Test that find_jpn returns correct values.
    """
    inputs = ["+10", [5, 4]], ["1.5-3", [9, 2]]
    expecteds = [1, 1, 10, False], [3, -1, 3, True]

    for i in range(len(inputs)):
        j, parity, n_states, is_jproj = split_jpn(*inputs[i])
    
        success = [j, parity, n_states, is_jproj] == expecteds[i]
        msg = "Bad return value from split_jpn."
        msg += f" Expected: {expecteds[i]}, got: {[j, parity, n_states, is_jproj]}."
        assert success, msg

if __name__ == "__main__":
    test_split_jpn_value_error()
    test_split_jpn()