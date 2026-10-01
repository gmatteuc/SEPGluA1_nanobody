function s = annotation_settings()
%ANNOTATION_SETTINGS  The settings of the control-point annotation.
%   s = ANNOTATION_SETTINGS() returns
%
%     s.suggested_anchors   how many anchor slices the GUI suggests for the
%                           automatic annotation, spread evenly from the
%                           first slice to the last (j jumps between them)
%     s.outlier_rule        when a pair counts as far off its slice's own
%                           affine, refitted without the outliers: a
%                           residual over .factor times the slice's median
%                           residual and over .min_px pixels; and when the
%                           slice as a whole is off: a median residual over
%                           .slice_px pixels (usually the wrong atlas plane)
%
%   One definition for the two places that use them: the control-point GUI,
%   through auto_annotation_plugin (the suggested anchors, the * mark on the
%   points and the slice title's warning), and register_to_atlas, which runs
%   the same check before registering (report_suspect_pairs). Without the
%   plugin, the GUI uses its own defaults, which are the values below.

s.suggested_anchors = 4;
s.outlier_rule = struct('factor', 3, 'min_px', 30, 'slice_px', 20);
end
