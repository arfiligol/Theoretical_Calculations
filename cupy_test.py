import cupy as cp

# 顯示 CUDA 驅動版本
print("CUDA 驅動版本:", cp.cuda.runtime.runtimeGetVersion())

# 顯示 CUDA 設備數量
device_count = cp.cuda.runtime.getDeviceCount()
print(f"可用的 CUDA 設備數量: {device_count}")

# 顯示 cuDNN 版本
print("cuDNN 版本:", cp.cuda.cudnn.getVersion())

# 顯示每個 CUDA 設備的信息
for i in range(device_count):
    device = cp.cuda.Device(i)
    print(f"設備 {i}: {device.name()}")
