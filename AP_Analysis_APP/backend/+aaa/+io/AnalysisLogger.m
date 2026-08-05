classdef AnalysisLogger < handle
    %ANALYSISLOGGER Timestamped plain-text logger for App and workflow use.

    properties (SetAccess = private)
        FilePath string = ""
        EchoToCommandWindow logical = true
    end

    properties (Access = private)
        FileId double = -1
    end

    methods
        function obj = AnalysisLogger(filePath, echoToCommandWindow)
            if nargin < 1
                filePath = "";
            end
            if nargin >= 2 && ~isempty(echoToCommandWindow)
                obj.EchoToCommandWindow = logical(echoToCommandWindow);
            end
            obj.FilePath = string(filePath);
            if strlength(obj.FilePath) > 0
                folder = fileparts(obj.FilePath);
                if ~isfolder(folder)
                    mkdir(folder);
                end
                obj.FileId = fopen(obj.FilePath, 'a', 'n', 'UTF-8');
                if obj.FileId < 0
                    error('AAA:IO:LogOpenFailed', ...
                        'Could not open analysis log: %s', obj.FilePath);
                end
            end
        end

        function write(obj, level, message)
            level = upper(strtrim(string(level)));
            message = string(message);
            if ~isscalar(message)
                message = strjoin(message(:)', newline);
            end
            line = sprintf('%s | %-7s | %s', ...
                char(datetime('now','Format','yyyy-MM-dd HH:mm:ss.SSS')), ...
                char(level), char(message));
            if obj.EchoToCommandWindow
                fprintf('%s\n', line);
            end
            if obj.FileId >= 0
                fprintf(obj.FileId, '%s\n', line);
                % MATLAB has no portable fflush. Close/reopen so App log
                % viewers can read each event immediately during long runs.
                fclose(obj.FileId);
                obj.FileId = fopen(obj.FilePath, 'a', 'n', 'UTF-8');
                if obj.FileId < 0
                    error('AAA:IO:LogReopenFailed', ...
                        'Could not reopen analysis log: %s', obj.FilePath);
                end
            end
        end

        function close(obj)
            if obj.FileId >= 0
                fclose(obj.FileId);
                obj.FileId = -1;
            end
        end

        function delete(obj)
            obj.close();
        end
    end
end
