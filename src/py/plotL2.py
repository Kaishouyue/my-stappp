import numpy as np
import matplotlib.pyplot as plt

# Element sizes and corresponding L2 error norms
h = [np.sqrt(2)/4, np.sqrt(2)/2, np.sqrt(2)]
L2 = [0.005059949324838033, 0.016013106482441183, 0.052956593577287416]

plt.figure(figsize=(6,4))
plt.loglog(h, L2, 'o-', label='L2 error')
plt.xlabel('Element size h')
plt.ylabel('L2 error norm')
plt.title('Log-Log Plot of L2 Error Norm vs. Element Size')
plt.grid(True, which="both", ls="--")
plt.legend()
plt.tight_layout()
plt.savefig("plotL2.png", dpi=300)
plt.show()

# 计算log-log斜率
log_h = np.log(h)
log_L2 = np.log(L2)
slope, intercept = np.polyfit(log_h, log_L2, 1)
print(f"Slope (convergence rate) = {slope}")
