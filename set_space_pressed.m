function set_space_pressed(event, fig)
    if strcmp(event.Key, 'space')
        fig.UserData = 'spacePressed';
    elseif strcmp(event.Key, 'return')
        fig.UserData = 'stop';
    end
end