function emit_phase(ctx,section,phase,state,message,role)
%EMIT_PHASE Emit a functional substage event without changing calculations.

if nargin<6,role="";end
event=struct('type',"phase",'mode',string(ctx.mode), ...
    'scope',string(ctx.scope),'input_path',string(ctx.input_path), ...
    'output_path',string(ctx.output_path),'section',string(section), ...
    'phase',string(phase),'role',string(role),'state',string(state), ...
    'message',string(message),'timestamp',datetime('now'));
if string(ctx.scope)=="single"
    [~,cycle_name]=fileparts(char(string(ctx.input_path)));
    if startsWith(string(cycle_name),"Cycle",'IgnoreCase',true)
        event.cycle=string(cycle_name);
    end
end
aaa.io.emit_event(ctx.runtime,event);
end
