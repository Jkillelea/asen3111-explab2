% lift slope = NACA 0012 slope
% No twist
% Note that I convert everything to feet from mm. I don't know if that's nessecary but
% when we wrote this code everything was in feet

% LONG
% chord = 60mm
% span = 240mm
% NACA 0012
long_chord = mm2ft(60);
long_span  = mm2ft(240);


% SHORT
% chord = 60mm
% span = 180mm
% NACA 0012
short_chord = mm2ft(60);
short_span  = mm2ft(180);

% ELIPTICAL
% root chord = 60mm
% tip chord = 1mm
% span = 180mm
% NACA 0012
e_span = mm2ft(180); 
e_c_r  = mm2ft(60);
e_c_t  = mm2ft(10);


lift_slope = 2*pi;
aero_twist = 0;
geo_twist  = 5; % AoA, degrees

[e, ~, ~] = PLLT(long_span, 2*pi, 2*pi, ...
                    long_chord, long_chord, ...
                    aero_twist, aero_twist, ...
                    geo_twist, geo_twist, 100);
fprintf('Long wing: e = %f\n', e);

[e, ~, ~] = PLLT(short_span, 2*pi, 2*pi, ...
                    short_chord, short_chord, ...
                    aero_twist, aero_twist, ...
                    geo_twist, geo_twist, 100);
fprintf('Short wing: e = %f\n', e);

[e, ~, ~] = PLLT(e_span, 2*pi, 2*pi, ...
                    e_c_t, e_c_r, ...
                    aero_twist, aero_twist, ...
                    geo_twist, geo_twist, 100);
fprintf('Elliptical wing: e = %f\n', e);
