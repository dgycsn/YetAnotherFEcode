clear
clc
close all

f1=300:3:310;
f2=310:0.5:320;
f3=320.25:0.25:325;
f4=327:2:335;
f5=340:5:350;

freq_forward=[f1 f2 f3 f4 f5];
freq_backward=flip(freq_forward);

force_amplitude=[1];
%%
% pause(10)
tic
[lin_forward_full, nl_forward_full,lin_forward_red,nl_forward_red,nl_forward_t]=...
    freq_resp_1d_new(freq_forward,force_amplitude);

[lin_backward_full, nl_backward_full,lin_backward_red,nl_backward_red,nl_backward_t]=...
    freq_resp_1d_new(freq_backward,force_amplitude);
total_time=toc;
%% plotting
figure;hold on

% plot(freq_forward,lin_forward_full,"DisplayName","lin forward full",...
%     "LineStyle","--","LineWidth",2,"Color","r","Marker","o")
% plot(freq_forward,nl_forward_full,"DisplayName","nonlin forward full",...
%     "LineStyle",":","LineWidth",2,"Color","r","Marker","s")
% plot(freq_forward,lin_forward_red,"DisplayName","lin forward red",...
%     "LineStyle","--","LineWidth",2,"Color","r","Marker","x")
% plot(freq_forward,nl_forward_red,"DisplayName","nonlin forward red",...
%     "LineStyle",":","LineWidth",2,"Color","r","Marker","d")
% plot(freq_forward,nl_forward_t,"DisplayName","nonlin forward t",...
%     "LineStyle","-.","LineWidth",2,"Color","r","Marker","*")
% 
% 
% plot(freq_backward,lin_backward_full,"DisplayName","lin backward full",...
%     "LineStyle","--","LineWidth",2,"Color","b","Marker","o")
% plot(freq_backward,nl_backward_full,"DisplayName","nonlin backward full",...
%     "LineStyle",":","LineWidth",2,"Color","b","Marker","s")
% plot(freq_backward,lin_backward_red,"DisplayName","lin backward red",...
%     "LineStyle","--","LineWidth",2,"Color","b","Marker","x")
% plot(freq_backward,nl_backward_red,"DisplayName","nonlin backward red",...
%     "LineStyle",":","LineWidth",2,"Color","b","Marker","d")
% plot(freq_backward,nl_backward_t,"DisplayName","nonlin backward t",...
%     "LineStyle","-.","LineWidth",2,"Color","b","Marker","*")

load("nlvib_freq_resp_30_8750.mat")
plot(nlvib_lin_freq,nlvib_lin_amp,"DisplayName","nlvib lin",...
    "LineStyle","--","LineWidth",2,"Color","k")
plot(nlvib_nonlin_freq,nlvib_nonlin_amp,"DisplayName","nlvib nonlin",...
    "LineStyle",":","LineWidth",2,"Color","k")
% f_s=freq_forward(1);
% f_e=freq_forward(end);
% xlim([f_s,f_e])
legend