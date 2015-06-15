/*
    David Jackson
*/
// selectableTable is an extension of jquery ui selectable that is faster on large tables
(function($) {

    $.widget("ui.selectableTable", $.ui.selectable, {
        _create: function() {
            var self = this;
            var first = true;
            var autoRefresh = this.options.autoRefresh;
	    bSort = false;
            this.element.addClass("ui-selectable");

            this.dragged = false;

            // cache selectee children based on filter
            this.refresh = function() {
                //console.log("ui.selectableTable.refresh()");
                var cellPositions = [];
                var firstRow = true;
                selecteeRows = $("tr", self.element[0]);
                selecteeRows.each(function() {
                    var row = $(this);
                    // adding any filters here seems to slow down IE a lot
                    // self.options.filter is not used!!
                    //var selecteeCells = $(row).find("td,th"); 
                    var selecteeCells = $(row).find("td.selectable, th.selectable");  // add filter here (jinlee)

                    // if (firstRow && selecteeCells.length > 0 ) {
                    if (firstRow && selecteeCells.length > 0 ) { 

                        for (var i = 0; i < selecteeCells.length; i++) {
                            cellPositions[i] = $(selecteeCells[i]).offset().left;
                        }
                        cellPositions[selecteeCells.length] = cellPositions[selecteeCells.length - 1] + $(selecteeCells).outerWidth(true);
                        firstRow = false;
                    }

                    if (selecteeCells.length > 0) {
                        var top = $(selecteeCells).first().offset().top;
                        var bottom = top + $(selecteeCells).first().outerHeight();
                        var i = 0;
                        selecteeCells.each(function() {
                            var $this = $(this);
                            first = false;
                            $.data(this, "selectable-item", {
                                element: this,
                                $element: $this,
                                left: cellPositions[i],
                                top: top,
                                right: cellPositions[i + 1],
                                bottom: bottom,
                                startselected: false,
                                selected: $this.hasClass('ui-selected'),
                                selecting: $this.hasClass('ui-selecting'),
                                unselecting: $this.hasClass('ui-unselecting')
                            });
                            i++;
                        });
                    }
                });

            };

            this.refresh();
            this.selectees = $(self.options.filter, self.element[0]).addClass("ui-selectee");

            this._mouseInit();

            this.helper = $(document.createElement('div'))
  		.css({ border: '1px dotted black' })
			.addClass("ui-selectable-helper");
        }
    });   
})(jQuery);
