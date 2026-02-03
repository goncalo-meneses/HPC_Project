all: build/power_gpu_32.x build/power_gpu_64.x build/power_gpu_128.x build/global_power_gpu.x

step_1: build/global_power_gpu.x build/power_gpu.x

step_2: build/power_gpu_32.x build/power_gpu_64.x build/power_gpu_128.x

build/power_gpu_32.x: power_gpu.cu | build
	nvcc -DBLOCK_SIZE=32 -o build/power_gpu_32.x power_gpu.cu

build/power_gpu_64.x: power_gpu.cu | build
	nvcc -DBLOCK_SIZE=64 -o build/power_gpu_64.x power_gpu.cu

build/power_gpu_128.x: power_gpu.cu | build
	nvcc -DBLOCK_SIZE=128 -o build/power_gpu_128.x power_gpu.cu

build/power_gpu.x: power_gpu.cu | build
	nvcc -DBLOCK_SIZE=32 -o build/power_gpu.x power_gpu.cu

build/global_power_gpu.x: global_power_gpu.cu | build
	nvcc -o build/global_power_gpu.x global_power_gpu.cu

build:
	mkdir -p build

clean:
	rm -rf build
