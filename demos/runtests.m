function runtests
% Run tests
tests = {'demo1', 'demo2', 'demo2a', 'demo3', 'demo3b', 'demo3c', 'demo4', 'demo5', 'demo6', ...
		'demo7', 'demo8', 'demo9', 'demo9a', 'demo11', 'geodemo_1a', 'geodemo_1b', 'geodemo_1c', ...
		'geodemo_1d', 'geodemo_1e', 'geodemo_2', 'geodemo_3','geodemo_6'};
		% broken URLs: demo10 geodemo_4 geodemo_4b
        % Not valid on master: access_test1
good = 0;
bad = 0;
failed = {};
for t = tests
	s = char(t);
	try
		fprintf(1, '\n\n==================================================================\n');
		fprintf(1, '=== %s ===========================================================\n', s);
		fprintf(1, '==================================================================\n');
		eval(s);
		close all;
		good = good + 1;
	catch me
		fprintf(1, '!!! %s failed: %s\n', s, me.identifier);
		fprintf(1, '%s\n', me.message);
		bad = bad + 1;
		failed{bad} = s;
	end
end
fprintf(1, '\n\n==========================================================\n');
fprintf(1, 'Ran %.0f demos\n', length(tests));
fprintf(1, '  %.0f passed\n', good);
fprintf(1, '  %.0f failed:\n', bad);
for f = failed
	fprintf(1, '      %s\n', char(f));
end
