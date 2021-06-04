function [out1,out2,out3]=create_settings_grid(in1,in2,in3)

out1=repelem(in1,numel(in2)*numel(in3));

out2_dm=repelem(in2,numel(in3));

out2=repmat(out2_dm,[1,numel(in1)]);

out3=repmat(in3,[1,numel(in1)*numel(in2)]);
