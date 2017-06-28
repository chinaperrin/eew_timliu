function gf = get_gain_factor_jp(scfString)
% Computes the factor that turns Japanes raw strong motion records [counts]
% into [ms^-2]
%
% --> s = sraw*gf;

slashIdx     = find(scfString=='/');
scaleFactor1 = str2double(scfString(slashIdx+1:end));
bracketIdx   = find(scfString=='(');
scaleFactor2 = str2double(scfString(1:bracketIdx-1));

gf = scaleFactor2/scaleFactor1*0.01;

