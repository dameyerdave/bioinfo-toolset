default:
	@echo "use deploy or clean"
clean:
	@rm -rf build dist
deploy: clean
	@python setup.py sdist bdist_wheel
	@twine upload dist/*