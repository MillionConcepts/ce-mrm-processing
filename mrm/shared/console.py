from sys import stdout


def print_inline(text: str, blanks: int = 60):
    """Classic utility function for quick and easy status displays."""
    stdout.write(" " * blanks + "\r")
    stdout.write(str(text) + "\r")
    stdout.flush()
