function h=redraw_text(h)

    props=get(h);
    delete(h);
    x=props.Position(1);
    y=props.Position(2);
    z=props.Position(3);

    h=text(x,y,z,'string','string',props.String,'color',props.Color);
    set(h,'BackgroundColor',props.BackgroundColor);
    set(h,'BusyAction',props.BusyAction);
    set(h,'ButtonDownFcn',props.ButtonDownFcn);
    set(h,'Children',props.Children);
    set(h,'Clipping',props.Clipping);
    set(h,'Color',props.Color);
    set(h,'CreateFcn',props.CreateFcn);
    set(h,'DeleteFcn',props.DeleteFcn);
    set(h,'DisplayName',props.DisplayName);
    set(h,'EdgeColor',props.EdgeColor);
    set(h,'Editing',props.Editing);
    set(h,'FontAngle',props.FontAngle);
    set(h,'FontName',props.FontName);
    set(h,'FontSize',props.FontSize);
    set(h,'FontUnits',props.FontUnits);
    set(h,'FontWeight',props.FontWeight);
    set(h,'HandleVisibility',props.HandleVisibility);
    set(h,'HitTest',props.HitTest);
    set(h,'HorizontalAlignment',props.HorizontalAlignment);
    set(h,'Interpreter',props.Interpreter);
    set(h,'Interruptible',props.Interruptible);
    set(h,'LineStyle',props.LineStyle);
    set(h,'LineWidth',props.LineWidth);
    set(h,'Margin',props.Margin);
    set(h,'Parent',props.Parent);
    set(h,'Position',props.Position);
    set(h,'Rotation',props.Rotation);
    set(h,'Selected',props.Selected);
    set(h,'SelectionHighlight',props.SelectionHighlight);
    set(h,'String',props.String);
    set(h,'Tag',props.Tag);
    set(h,'UIContextMenu',props.UIContextMenu);
    set(h,'Units',props.Units);
    set(h,'UserData',props.UserData);
    set(h,'VerticalAlignment',props.VerticalAlignment);
    set(h,'Visible',props.Visible);

