import numpy

def get_dielectic(lamb):
    with open('refractiveindex/esAG1.txt','r+') as file:
        lambs = list(map(lambda x: float(x), file.readline().split(',')))
        epsilons = list(map(lambda x: complex(x), file.readline().replace('i', 'j').split(',')))
        epsilon=epsilons[0]
        if len(lambs) == len(epsilons):
            for i in range(0, len(lambs)):
                if lambs[i-1] <= lamb*(10**(-6)) <= lambs[i]:
                    if (lambs[i-1] + abs((lambs[i] - lambs[i-1]) / 2)) >= lamb:
                        epsilon = epsilons[i-1]
                    else:
                        epsilon = epsilons[i]

    return epsilon

print(get_dielectic(0.1879))

#-0.32404+2.5937i
#1.879e-07

#-0.32324+2.597i
#1.88e-07
