function runtests
% Run tests
good = 0;
bad = 0;
failed = {};
mydir= dir('*geodemo*.m');
for itest=1:length(mydir);
    ss = char(mydir(itest).name(1:end-2))
    tests{itest}=ss;
	try
		fprintf(1, '\n\n===============================================================\n')
		fprintf(1, '=== %s =====================================================\n', ss)
		fprintf(1, '===============================================================\n')
		eval(ss);
		close all;
		good = good + 1;
	catch me
		fprintf(1, '!!! %s failed: %s\n', ss, me.identifier);
		fprintf(1, '%s\n', me.message);
		bad = bad + 1;
		failed{bad} = ss;
	end
end
fprintf(1, '\n\n==========================================================\n')
fprintf(1, 'Ran %.0f demos\n', length(tests))
fprintf(1, '  %.0f passed\n', good)
fprintf(1, '  %.0f failed:\n', bad)
for f = failed
	fprintf(1, '      %s\n', char(f))
end
