{% - Extends 'slides_reveal.tpl' -%}

{% Block input_group -%}
<Div class = "input_hidden">
{{Super ()}}
</ Div>
{% Endblock%} input_group

{% - Block header -%}
{{Super ()}}

<Script src = "// ajax.googleapis.com/ajax/libs/jquery/1.10.2/jquery.min.js"> </ script>

<Style type = "text / css">
//div.output_wrapper {
// Margin-top: 0px;
} //
.input_hidden {
  display: none;
// Margin-top: 5px;
}
</ Style>

<Script>
$ (Document) .ready (function () {
  $ (". Output_wrapper"). Click (function () {
      $ (This) .prev ('input_hidden.') SlideToggle ().;
  });
})
</ Script>
{% - Endblock header -%}
