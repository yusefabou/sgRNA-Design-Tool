from django.conf.urls import url

from . import views

urlpatterns = [
	url(r'^$', views.index, name='index'),
	url(r'^form.html$', views.form, name='form'),
	url(r'^result.html$', views.result, name='result')
]