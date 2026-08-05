function files=render_registration_preview(preview,mode,output_path)
%RENDER_REGISTRATION_PREVIEW Save the frozen mode-specific registration QC.
files=struct();mode=lower(string(mode));if mode=="none"||isempty(fieldnames(preview)),return;end
fig=figure('Color','w','Visible','off');subplot(1,3,1);imagesc(preview.voltage_mean);axis image off;title('Voltage average');subplot(1,3,2);imagesc(preview.calcium_mean);axis image off;title('Calcium average');subplot(1,3,3);imshowpair(preview.voltage_mean,preview.calcium_registered);title('Registered overlay');
if mode=="manual_points",stem="0_dual_manual_points_preview";else,stem="0_dual_matlab_register_preview";end
files=struct('fig_file',string(fullfile(output_path,stem+".fig")),'png_file',string(fullfile(output_path,stem+".png")));savefig(fig,files.fig_file);exportgraphics(fig,files.png_file,'Resolution',150);close(fig);
end
