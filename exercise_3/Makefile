all: build/power_gpu.x build/global_power_gpu.x

build/power_gpu.x: power_gpu.cu | build
	nvcc -o build/power_gpu.x power_gpu.cu

build/global_power_gpu.x: global_power_gpu.cu | build
	nvcc -o build/global_power_gpu.x global_power_gpu.cu

build:
	mkdir -p build

clean:
	rm -rf build
