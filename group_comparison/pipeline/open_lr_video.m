function [vidObj, full_video_path] = open_lr_video(save_dir, video_filename)
%OPEN_LR_VIDEO  Open the MPEG-4 file of a left-right video, its folder made if needed.
%   [vidObj, full_video_path] = OPEN_LR_VIDEO(save_dir, video_filename) makes
%   save_dir if it is missing, and opens save_dir\video_filename for writing at
%   15 frames per second, quality 95.
%
%   Run by the write_lr_* video writers.

% create the folder if needed
if ~exist(save_dir, 'dir')
    mkdir(save_dir);
end

% open the video
full_video_path = fullfile(save_dir, video_filename);
vidObj = VideoWriter(full_video_path, 'MPEG-4');
vidObj.FrameRate = 15;
vidObj.Quality = 95;
open(vidObj);

end
