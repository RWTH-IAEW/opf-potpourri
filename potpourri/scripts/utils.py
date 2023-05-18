import os
import platform
import socket
from datetime import datetime
import shutil

# ANSI escape sequences for color formatting
color_red = "\033[91m"
color_green = "\033[92m"
color_yellow = "\033[93m"
color_reset = "\033[0m"

def endLine(color: str = None) -> None:
    import shutil

    # Get the width of the terminal console
    terminal_width = shutil.get_terminal_size().columns

    # Define the decorative bar pattern
    bar_pattern = "\u2500"  # Customize the pattern as desired

    # Create the decorative horizontal bar
    horizontal_bar = bar_pattern * terminal_width
    
    if color == "red":
        print(color_red + horizontal_bar + color_reset)
    elif color == "green":
        print(color_green + horizontal_bar + color_reset)
    elif color == "yellow":
        print(color_yellow + horizontal_bar + color_reset)
    else:
        print(horizontal_bar)
    
    return

def newHeading(heading_text:str = "New Section") -> None:

    # Calculate the total width of the heading
    heading_width = len(heading_text) + 4  # Add extra space for the borders

    # Create the Unicode box characters
    horizontal_line = "\u2500" * heading_width
    vertical_line = "\u2502"
    top_border = "\u250C" + horizontal_line + "\u2510"
    bottom_border = "\u2514" + horizontal_line + "\u2518"

    # Print the heading
    print(top_border)
    print(vertical_line, heading_text, vertical_line)
    print(bottom_border)
    
    return

def showParameter(testcase:str) -> None:
    endLine()
    print(color_green + " ____       _                               _ " + color_reset)
    print(color_green + "|  _ \ ___ | |_ _ __   ___  _   _ _ __ _ __(_)" + color_reset)
    print(color_green + "| |_) / _ \| __| '_ \ / _ \| | | | '__| '__| |" + color_reset)
    print(color_green + "|  __/ (_) | |_| |_) | (_) | |_| | |  | |  | |" + color_reset)
    print(color_green + "|_|   \___/ \__| .__/ \___/ \__,_|_|  |_|  |_|" + color_reset)
    print(color_green + "               |_|                            " + color_reset)
    print("TESTCASE:  " + testcase) #rural_test
    print("RUN:       " + datetime.now().strftime("%d/%m/%Y %H:%M:%S")) #29 - Apr - 2022    08: 54:37
    print("OS:        " + platform.system())
    print("OM:        " + socket.gethostname())  # PC4195
    if platform.system() not in ["linux", "linux2", "Linux", "Linux2"]:
        print("USER:  " + os.environ.get("USERNAME"))# pr092066
    endLine()
    
    return

if __name__ == "__main__":
    # Example output with colored text
    print(color_red + "This text is red." + color_reset)
    print(color_green + "This text is green." + color_reset)
    print(color_yellow + "This text is yellow." + color_reset)

