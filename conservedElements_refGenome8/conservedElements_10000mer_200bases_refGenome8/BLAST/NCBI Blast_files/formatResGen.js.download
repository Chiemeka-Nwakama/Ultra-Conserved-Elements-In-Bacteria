// JScript source code
//takes care of jQuery calls to elements that have dots in id like id= "SRA34456.3.4 or SRA:ERR1293056.5.1"
function initjQry() 
{
    var oldCreate = jQuery.ui.ncbipopper.prototype._create;
    jQuery.ui.ncbipopper.prototype._create = function () {
        var destSelector = this.options.destSelector || this.options.sourceSelector || this.element.prop("hash");
        if (typeof destSelector == "string" ) {
            destSelector = destSelector.replace(/\./g, "\\.");
            destSelector = destSelector.replace(/\:/g, "\\:");
            this.options.destSelector = destSelector;
            //this.options.destSelector = destSelector.replace(/\./g, "\\.");
            oldCreate.apply(this, arguments);
        }
    }    
}

function initResultspage(e)
{
    initAdvancedView(e);    
    if ($("smrtBlastFullSearch")) initSmartBlast();
    if ($("blastpSubmit")) jQuery("#blastpSubmit").bind("click", function (e) { submitCurrOptions(this); });
    jQuery("#grBlastHits").find("a.graphSeq").each(function (index) {        
        jQuery(this).bind("ncbipopperopen", function (e) {
            var seq = this.getAttribute("seqid");
            var score = this.getAttribute("score");
            var eval = this.getAttribute("eval");
            jQuery("#dlgGraph_" + seq).find("span.sc").each(function (index) {
                jQuery(this).html(score);
            });
            jQuery("#dlgGraph_" + seq).find("span.ev").each(function (index) {
                jQuery(this).html(eval);
            });
        });
    });
    jQuery("#changeView").bind("click", function(e) {               
       var paramString = "ADV_VIEW=" + jQuery(this).attr("advview") + "&CONFIG_DESCR=" + jQuery(this).attr("configdescr");       
       this.href += "&" + paramString;       
       SaveOptionsCubby(paramString);
    });        
    displayResultsPageSurvey("tRaiETqnLgj758hTBazgdwxw_2BrMG_2BYwhTjHo2Ki3Pz5Eieb_2BAj1z5J9jsyPzOmi2");
}
				
utils.addEvent(window,"load", initResultspage,false);   

initjQry();
