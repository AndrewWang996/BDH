function r = hasGPUComputing()

persistent hasgpu;
if isempty(hasgpu), hasgpu = gpuDeviceCount>0; end

r = hasgpu;