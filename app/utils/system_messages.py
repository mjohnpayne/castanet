from termcolor import colored


def banner():
    print("""
 _____           _                   _
/  __ \         | |                 | |
| /  \/ __ _ ___| |_ __ _ _ __   ___| |_
| |    / _` / __| __/ _` | '_ \ / _ \ __|
| \__/\ (_| \__ \ || (_| | | | |  __/ |_
 \____/\__,_|___/\__\__,_|_| |_|\___|\__|
O       o O       o O       o O       o O
| O   o | | O   o | | O   o | | O   o | |
| | O | | | | O | | | | O | | | | O | | |
| o   O | | o   O | | o   O | | o   O | |
o       O o       O o       O O       O o
(c) Rich Mayne & Tanya Golubchik
    https://doi.org/10.1093/bioinformatics/btae591
""")


def end_sec_print(msg):
    print(f"\n{'*'*30}\n{colored(msg, 'green')}\n{'*'*30}\n")
