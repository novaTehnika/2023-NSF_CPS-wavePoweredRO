%% into RO
val = 1e3*0.0112*(out.q_feed - out.q_ERUfeed);
min(val)
mean(val)
max(val)

%% into RO accumulators
val = 1e3*0.0112*abs(max(0,out.q_rv - (out.q_feed - out.q_ERUfeed)));
min(val)
mean(val)
max(val)
prctile(val,97)

% out of RO accumulators
val = 1e3*0.0112*abs(min(0,out.q_rv - (out.q_feed - out.q_ERUfeed)));
min(val)
mean(val)
max(val)
prctile(val,97)

%% into WEC-driven Pump accumulators
val = 1e3*0.0112*abs(max(0,out.q_h - out.q_pm - out.q_rv));
min(val)
mean(val)
max(val)
prctile(val,97)

% out of WEC-driven Pump accumulators
val = 1e3*0.0112*abs(min(0,out.q_h - out.q_pm - out.q_rv));
min(val)
mean(val)
max(val)
prctile(val,97)

%% into LP accumulators
val = 1e3*0.0112*abs(max(0,(out.q_c - out.q_ERUfeed) - out.q_h + out.q_pm));
min(val)
mean(val)
max(val)

% out of LP accumulators
val = 1e3*0.0112*abs(min(0,(out.q_c - out.q_ERUfeed) - out.q_h + out.q_pm));
min(val)
mean(val)
max(val)

%% charge Pump
val = 1e3*0.0112*(out.q_c - out.q_ERUfeed);
min(val)
mean(val)
max(val)

%% WEC-driven Pump
val = 1e3*0.0112*out.q_h;
min(val)
mean(val)
max(val)


%% Hyd. motor/gen
val = 1e3*0.0112*out.q_pm;
min(val)
mean(val)
max(val)

%% WEC-driven pump minus hyd. motor
val = 1e3*0.0112*(out.q_h - out.q_pm);
min(val)
mean(val)
max(val)


%%
q = 0.49/1e3;
d = 0.375*2.54/100;
A = pi*d^2/4;
v_mean = q/A