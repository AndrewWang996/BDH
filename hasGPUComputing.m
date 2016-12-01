function r = hasGPUComputing()
fprintf('\n Running hasGPUComputing');
persistent hasgpu;
if isempty(hasgpu), hasgpu = gpuDeviceCount>0; end

r = hasgpu;