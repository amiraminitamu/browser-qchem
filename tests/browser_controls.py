"""Reveal a workbench popover or drawer before a browser test edits its control."""
def control(page, selector):
    locator = page.locator(selector)
    parents = locator.evaluate('''element => {
        const panels=[];
        for(let parent=element.parentElement;parent;parent=parent.parentElement)
            if(parent.tagName==='DETAILS' && parent.id) panels.unshift(parent.id);
        return panels;
    }''')
    for id in parents:
        panel = page.locator('#'+id)
        if not panel.evaluate('(element)=>element.open'):
            if id == 'coordinates-editor':
                page.locator('#xyz-toggle').click()
            else:
                panel.locator(':scope > summary').click()
    return locator
