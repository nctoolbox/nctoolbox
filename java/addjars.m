function addjars(jarDir)

javaaddpath(jarDir);

d = dir(jarDir);
for i = 1:length(d)
    name = java.lang.String(d(i).name);
    if name.endsWith('.jar') & ~name.startsWith('.')
        javaaddpath(fullfile(jarDir, d(i).name));
    end
end