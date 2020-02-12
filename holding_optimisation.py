import sys


def readInput(data):
    c1, c2, sigma_h, h_a_0, h_b_0 = data[0].split(" ")
    return c1, c2, sigma_h, h_a_0, h_b_0

def main():
    data = sys.stdin.readlines()
    c1, c2, sigma_h, h_a_0, h_b_0 = readInput(data)
    print (h_b_0)

main()

