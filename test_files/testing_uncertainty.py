import matplotlib.pyplot as plt
import numpy as np
import math

from matplotlib import animation


def prob_func_phi_est(fA,fB,A,B,C,phi_A,phi_B,phi_C):
    #fA = px  # check_point[0] / (check_point[0] + check_point[1])  # sum(check_point)
    #fB = py  # check_point[1] / (check_point[0] + check_point[1])  # sum(check_point)

    x = phi_A
    y = phi_B
    z = phi_C

    measured_point = (A, B, C)
    print(measured_point)
    mfA = measured_point[0] / sum(measured_point)
    mfB = measured_point[1] / sum(measured_point)
    sigfA = math.sqrt((x**2*(B*C*(y + z)**2 + A*(C*y**2 + B*z**2)))/(A*B*C*(x + y + z)**4))  # math.sqrt(mfA * (1 - mfA) / sum(measured_point))
    sigfB = math.sqrt((y**2*(B*C*x**2 + A*(B*z**2 + C*(x + z)**2)))/(A*B*C*(x + y + z)**4))  # math.sqrt(mfB * (1 - mfB) / sum(measured_point))
    sigAB = -((x*y*(-(A*B*z**2) + A*C*y*(x + z) + B*C*x*(y + z)))/(A*B*C*(x + y + z)**4))
    rho = sigAB / (sigfA * sigfB)
    p1 = 1 / (2 * math.pi * sigfA * sigfB * np.sqrt(1 - rho * rho))
    p3 = (fA - mfA) * (fA - mfA) / (sigfA * sigfA) + (fB - mfB) * (fB - mfB) / (sigfB * sigfB) - 2 * rho * (
            fA - mfA) * (fB - mfB) / (sigfA * sigfB)

    p2 = np.exp(-.5 / (1 - rho * rho) * p3)
    prob = p1 * p2
    prob[prob > 1] = 1
    return prob

def prob_func(px,py,A,B,C):
    # global A
    # global B
    # global C
    fA = px # check_point[0] / (check_point[0] + check_point[1])  # sum(check_point)
    fB = py #check_point[1] / (check_point[0] + check_point[1])  # sum(check_point)
    measured_point = (A,B,C)
    print(measured_point)
    mfA = measured_point[0] / sum(measured_point)
    mfB = measured_point[1] / sum(measured_point)
    sigfA = math.sqrt(A*(B+C)/((A+B+C)**3)) #math.sqrt(mfA * (1 - mfA) / sum(measured_point))
    sigfB = math.sqrt(B*(A+C)/((A+B+C)**3)) #math.sqrt(mfB * (1 - mfB) / sum(measured_point))
    sigAB = -1 * A * B / math.pow(sum(measured_point), 3)
    rho = sigAB / (sigfA * sigfB)
    p1 = 1 / (2 * math.pi * sigfA * sigfB * np.sqrt(1 - rho * rho))
    p3 = (fA - mfA) * (fA - mfA) / (sigfA * sigfA) + (fB - mfB) * (fB - mfB) / (sigfB * sigfB) - 2 * rho * (
                fA - mfA) * (fB - mfB) / (sigfA * sigfB)

    p2 = np.exp(-.5 / (1 - rho * rho) * p3)
    prob = p1 * p2
    prob[prob > 1] = 1

    exp_expression = -1/(2*(1-rho**2))*(((fA-mfA)**2 / sigfA**2) + ((fB-mfB)**2 / sigfB**2) - (2*rho*(fA-mfA)*(fB-mfB)/(sigfA*sigfB)))
    prob2 = 1/(2*math.pi*sigfA*sigfB)*(1/math.sqrt(1-rho**2))*np.exp(exp_expression)
    prob2[prob2>1] = 1

    # return px**2 + py**2
    prob1 = ((A + B + C)**5*np.exp(-0.5*((A + B + C)*(B*C*(B + C)*px**2 + A**2*(C*py**2 + B*(-1 + px + py)**2) +\
                A*(C**2*py**2 + B**2*(-1 + px + py)**2 + 2*B*C*(px**2 + px*(-1 + py) + (-1 + py)*py))))/(A*B*C)))/(A*B*C*math.pi)
    prob1[prob1 > 1] = 1
    return prob2

if __name__ == '__main__':
    # see if we can just plot some poisson distributions for A, B, and C
    # samples = 1000
    # A_poisson_samples = np.random.poisson(A, samples)
    # B_poisson_samples = np.random.poisson(B, samples)
    # C_poisson_samples = np.random.poisson(C, samples)
    #
    # poissons_fig, pois_ax = plt.subplots(1,1)
    # pois_ax.hist(A_poisson_samples, label='A')
    # pois_ax.hist(B_poisson_samples, label='B')
    # pois_ax.hist(C_poisson_samples, label='C')
    # poissons_fig.legend()
    # plt.show()

    # now, plot the PDF, though this will be in 3D, so use a mesh
    pdf_fig, pdf_ax = plt.subplots(1,1) # subplot_kw={'projection':'3d'}
    x = np.linspace(0,1,1000)
    y = np.linspace(0,1,1000)
    px, py = np.meshgrid(x,y)
    # pdf_ax.set_zlim(0,1)
    A = 1
    B = 1
    C = 1
    x = y = z = 1
    surf = pdf_ax.pcolormesh(px, py, np.array(prob_func_phi_est(px, py,A,B,C,x,y,z)))
    pdf_fig.colorbar(surf)
    pdf_ax.set_xlabel('fA')
    pdf_ax.set_ylabel('fB')

    # let's do an animation to show sensitivity to constants
    constants = np.linspace(2,50,48)
    def __update(const):
        surf = pdf_ax.pcolormesh(px, py, np.array(prob_func_phi_est(px, py, const, const, const, const, const, const)))
        pdf_ax.set_title(f'const={const}')
        # pdf_fig.colorbar(surf)


    ani = animation.FuncAnimation(pdf_fig, __update, constants, repeat=True)
    # code courtesy of Josh Q.
    print('Saving spectrogram video animation...')
    writermp4 = animation.FFMpegWriter(fps=5)
    ani.save('video.mp4', writer=writermp4)
    print('...Done')
    # pdf_fig.show()
