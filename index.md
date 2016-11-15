---
layout: default
title: Software
---

### trajectory alignment software

Here, you can [Download](https://github.com/apicco/trajectory_alignment/archive/master.zip) the software from the [github repository](https://github.com/apicco/trajectory_alignment/).

Read the [domcumentation](wiki/Home)

<ul class="wikiMenu">
	{% for p in site.pages %}
		{% if p.menu == "wiki" %}
			{% if p.wikiPageName == "Home" %}
	        <li><a class="post-link" href="{{ p.url | prepend: site.baseurl }}"> Documentation </a></li>
	        {% endif %}
	    {% endif %}
	{% endfor %}
</ul>.

References:

[Picco, _et al._, eLife 2015](http://dx.doi.org/10.7554/eLife.04535)

[Kukulski, Picco, _et al._, eLife 2016](http://dx.doi.org/10.7554/eLife.16036)
